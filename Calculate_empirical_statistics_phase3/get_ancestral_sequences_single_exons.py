#This is to calcualte divergence per site in exons and in their nearby intergenic regions. 
#This script directly uses the fasta file of human ancestor of GRCH37

import sys

s_chrom = sys.argv[1]                                                             
                                                                                
#get chromosome names                                                           
d_chrom_names = {}                                                              
f_info = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_info.txt", 'r')
for line in f_info:                                                             
    line1 = line.strip('\n')                                                    
    line2 = line1.split('\t')                                                   
    d_chrom_names[line2[6]] = line2[9]                                          
f_info.close()

#get coordinates of exons:                                                      
f_exons = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon2-6kb_CNE500_rec_numbp50_final_div.tab", 'r')
d_col = {}                                                                      
l_exons = []                                                                    
exon_num = 1                                                                    
d_chrom, d_strand, d_numbp50 = {}, {}, {}

d_exon_start, d_exon_end = {}, {}
d_5p_start, d_5p_end, d_3p_start, d_3p_end = {}, {}, {}, {}
for line in f_exons:                                                            
    lineA = line.strip('\n')                                                    
    lineB = lineA.split('\t')                                                   
    if lineB[0] == "exon":                                                      
        col_num = 0                                                             
        for x in lineB:                                                         
            d_col[x] = col_num                                                  
            col_num += 1                                                        
    else:                                                                       
        s_gene = lineB[d_col["gene"]]
        s_exon = "exon" + str(exon_num) + "_" + s_gene
        d_chrom[s_exon] = d_chrom_names[lineB[d_col["chrom"]]]   
        if d_chrom[s_exon] == s_chrom:
            l_exons.append(s_exon)
            d_strand[s_exon] = lineB[d_col["strand"]]
            d_exon_start[s_exon] = lineB[d_col["exon_start"]]
            d_exon_end[s_exon] = lineB[d_col["exon_end"]]
            d_numbp50[s_exon] = lineB[d_col["numbp50"]]
            if d_strand[s_exon] == "+":
                d_5p_end[s_exon] = int(d_exon_start[s_exon]) - 1
                d_5p_start[s_exon] = int(d_exon_start[s_exon]) - 4*int(d_numbp50[s_exon])
                d_3p_start[s_exon] = int(d_exon_end[s_exon]) + 1
                d_3p_end[s_exon] = int(d_exon_end[s_exon]) + 4*int(d_numbp50[s_exon])
            elif d_strand[s_exon] == "-":
                d_5p_start[s_exon] = int(d_exon_end[s_exon]) + 1
                d_5p_end[s_exon] = int(d_exon_end[s_exon]) + 4*int(d_numbp50[s_exon])
                d_3p_end[s_exon] = int(d_exon_start[s_exon]) - 1
                d_3p_start[s_exon] = int(d_exon_start[s_exon]) - 4*int(d_numbp50[s_exon])
        exon_num += 1                                                           
f_exons.close()
#print(d_chrom)
d_anc_exon, d_anc_5p, d_anc_3p = {}, {}, {}
for s_exon in l_exons:
    d_anc_exon[s_exon] = ""
    d_anc_5p[s_exon] = ""
    d_anc_3p[s_exon] = ""

#read in the ancestral sequences:
f_fa = open("/scratch/pjohri1/BgsDfeDemo_Human/HumanData/divergence/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_" + s_chrom.replace("chr", "") + ".fa", 'r')
line_start_coord = 0
for line in f_fa:
    line1 = line.strip('\n')
    if line1[0] != ">":
        if line_start_coord == 0:
            line_start_coord = 1
        else:
            line_start_coord += 100
        line_end_coord = line_start_coord + 100
        for s_exon in l_exons:
            if line_start_coord in range(int(d_exon_start[s_exon]), int(d_exon_end[s_exon])+1) or line_end_coord in range(int(d_exon_start[s_exon]), int(d_exon_end[s_exon])+1):#exon lies here:
                #start recording sequence:
                s_coord = line_start_coord
                for x in line1:
                    if s_coord in range(int(d_exon_start[s_exon]), int(d_exon_end[s_exon])+1):
                        d_anc_exon[s_exon] = d_anc_exon[s_exon] + x
                    s_coord += 1
            #5prime:
            if line_start_coord in range(int(d_5p_start[s_exon]), int(d_5p_end[s_exon])+1) or line_end_coord in range(int(d_5p_start[s_exon]), int(d_5p_end[s_exon])+1):#5' intergenic lies here
                #start recording sequence:
                s_coord = line_start_coord 
                for x in line1:
                    if s_coord in range(int(d_5p_start[s_exon]), int(d_5p_end[s_exon])+1):
                        d_anc_5p[s_exon] = d_anc_5p[s_exon] + x
                    s_coord += 1
            #3prime:
            if line_start_coord in range(int(d_3p_start[s_exon]), int(d_3p_end[s_exon])+1) or line_end_coord in range(int(d_3p_start[s_exon]), int(d_3p_end[s_exon])+1):#3' intergenic lies here
                #start recording sequence:
                s_coord = line_start_coord
                for x in line1:
                    if s_coord in range(int(d_3p_start[s_exon]), int(d_3p_end[s_exon])+1):
                        d_anc_3p[s_exon] = d_anc_3p[s_exon] + x
                    s_coord += 1
f_fa.close()

#write the results:
result_exon = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/divergence/ANC/exon/" + s_chrom + ".fa", 'w+')
result_5p = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/divergence/ANC/5p/" + s_chrom + ".fa", 'w+')
result_3p = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/divergence/ANC/3p/" + s_chrom + ".fa", 'w+')
for s_exon in l_exons:
    result_exon.write(">" + s_exon +"_exon_anc" + '\n')
    result_exon.write(d_anc_exon[s_exon] + '\n')
    #writing intergenic sequences with noncoding and then exonic
    if d_strand[s_exon] == "+":
        result_5p.write(">" + s_exon +"_5p_anc" + '\n')
        result_5p.write(d_anc_5p[s_exon] + '\n') #written in forward
        result_3p.write(">" + s_exon +"_3p_anc" + '\n')
        result_3p.write(d_anc_3p[s_exon][::-1] + '\n') #written in reverse
    elif d_strand[s_exon] == "-":
        result_5p.write(">" + s_exon +"_5p_anc" + '\n')
        result_5p.write(d_anc_5p[s_exon][::-1] + '\n') #written in reverse
        result_3p.write(">" + s_exon +"_3p_anc" + '\n')
        result_3p.write(d_anc_3p[s_exon] + '\n') #written in forward
result.close()

print ("done")

