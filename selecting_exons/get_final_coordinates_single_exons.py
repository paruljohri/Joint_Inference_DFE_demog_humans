#This is just to write coordinates of all exons once:

import sys

#get chromosome names                                                           
d_chrom_num = {}                                                                 
f_info = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_info.txt", 'r')
for line in f_info:                                                             
    line1 = line.strip('\n')                                                    
    line2 = line1.split('\t')                                                   
    d_chrom_num[line2[6]] = line2[9]                                             
f_info.close()

#get the coordinates of the exons:                                              
print("reading exon coordinates")                                               
f_exons = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon2-6kb_CNE500_rec_numbp50_final_div.tab", 'r')
d_col = {}                                                                      
l_exons = []                                                                    
exon_num = 1                                                                    
d_chr = {}                                                                      
d_strand, d_numbp50 = {}, {}                                                    
#d_exon_start, d_exon_end = {}, {}                                              
d_coord_start = {}                                                              
d_coord_end = {}                                                                
l_regions = ["exon", "5p", "3p"]                                                
for region in l_regions:                                                        
    d_coord_start[region] = {}                                                  
    d_coord_end[region] = {}
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
        s_chr = d_chrom_num[lineB[d_col["chrom"]]]                              
                                                           
        l_exons.append(s_exon)                                              
        d_chr[s_exon] = s_chr                                               
        d_strand[s_exon] = lineB[d_col["strand"]]                            
        d_coord_start["exon"][s_exon] = int(lineB[d_col["exon_start"]])     
        d_coord_end["exon"][s_exon] = int(lineB[d_col["exon_end"]])         
        d_numbp50[s_exon] = int(lineB[d_col["numbp50"]])                    
        exon_num += 1                                                           
f_exons.close()

#Get coordinates of different regions:                                          
d_direction = {}
for region in l_regions:
    d_direction[region] = {}
for exon in l_exons:
    d_direction["exon"][exon] = "forward"
for exon in l_exons:                                                            
    if d_strand[exon] == "+":                                                   
        d_coord_start["5p"][exon] = int(d_coord_start["exon"][exon]) - 4.0*int(d_numbp50[exon])
        d_coord_end["5p"][exon] = int(d_coord_start["exon"][exon]) - 1          
        d_direction["5p"][exon] = "forward"
        d_coord_start["3p"][exon] = int(d_coord_end["exon"][exon]) + 1          
        d_coord_end["3p"][exon] = int(d_coord_end["exon"][exon]) + 4.0*int(d_numbp50[exon])
        d_direction["3p"][exon] = "reverse"
    elif d_strand[exon] == "-":                                                 
        d_coord_start["5p"][exon] = int(d_coord_end["exon"][exon]) + 1          
        d_coord_end["5p"][exon] = int(d_coord_end["exon"][exon]) + 4.0*int(d_numbp50[exon])
        d_direction["5p"][exon] = "reverse"
        d_coord_start["3p"][exon] = int(d_coord_start["exon"][exon]) - 4.0*int(d_numbp50[exon])
        d_coord_end["3p"][exon] = int(d_coord_start["exon"][exon]) - 1          
        d_direction["3p"][exon] = "forward"
#write it all out:
result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/single_exon_coordinates.txt", 'w+')
result.write("exon" + '\t' + "region" + '\t' + "chrom" + '\t' + "start" + '\t' + "end" + '\t' + "strand" + '\t' + "direction" + '\n')
for exon in l_exons:
    for region in l_regions:
        result.write(exon + '\t' + region + '\t' + d_chr[exon] + '\t' + str(round(d_coord_start[region][exon])) + '\t' + str(round(d_coord_end[region][exon])) + '\t' + d_strand[exon] + '\t' + d_direction[region][exon] + '\n')

result.close()
print("done")


