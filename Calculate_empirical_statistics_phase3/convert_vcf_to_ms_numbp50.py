#This is to convert the vcf files into ms files:
#Need to align your ms such that it's intergenic first and then the exon:

import sys
s_yri_folder = sys.argv[1] #YRI/YRI_50
subfolder="Numbp50"

#get the strands and list of each exon:
f_exons = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon2-6kb_CNE500_rec_numbp50_final_div.tab", 'r')
d_col = {}                                                                      
l_exons = []                                                                    
exon_num = 1
d_strand, d_numbp50 = {}, {}
d_exon_start, d_exon_end = {}, {}
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
        l_exons.append(s_exon)
        d_strand[s_exon] = lineB[d_col["strand"]]
        d_exon_start[s_exon] = lineB[d_col["exon_start"]]
        d_exon_end[s_exon] = lineB[d_col["exon_end"]]
        d_numbp50[s_exon] = lineB[d_col["numbp50"]]
        exon_num += 1
f_exons.close()
#print (l_exons)

def convert_vcf_to_ms(f_VCF):
    l_posns = []
    d_haps = {}
    for line in f_VCF:
        line1 = line.strip('\n')
        if line1[0] == "#":
            line2 = line1.split('\t')
            indv_num = 1
            for x in line2[9:]:
                d_haps[indv_num] = ""
                d_haps[indv_num+1] = ""
                indv_num += 2
        if line1[0]!= "#":
            if "VT=SNP" in line1:#represents a SNP
                line2 = line1.split('\t')
                #check if the site is monomorphic:
                s_gts = ''.join(line2[9:])
                if "1" in s_gts:#site is polymorphic in YRI
                    l_posns.append(int(line2[1]))
                    indv_num = 1
                    for x in line2[9:]:
                        x1 = x.split("|")
                        d_haps[indv_num] = d_haps[indv_num] + x1[0]
                        d_haps[indv_num+1] = d_haps[indv_num+1] + x1[1]
                        indv_num += 2
    return(d_haps, l_posns)

def normalize_positions(l_positions, posn_start, posn_end, direction):
    l_positions_norm = []
    LENGTH = int(posn_end) - int(posn_start) + 1
    if direction == "forward":
        for posn in l_positions:
            l_positions_norm.append((int(posn)-int(posn_start))/float(LENGTH-1))#to match SLiM's numbering
    elif direction == "reverse":
        for posn in l_positions:
            l_positions_norm.append((int(posn_end)-int(posn))/float(LENGTH-1))#to match SLiM's numbering
    return(l_positions_norm)

def write_ms_file(result, l_posn_norm, d_hap_exon, direction):
    result.write("//" + '\n')
    result.write("segsites: " + str(len(l_posn_norm)) + '\n')
    if direction == "forward":
        result.write("positions: " + ' '.join(str(e) for e in l_posn_norm) + '\n')
        for indv_num in range(1,len(d_hap_exon)+1):
            result.write(d_hap_exon[indv_num] + '\n')
    elif direction == "reverse":
        l_posn_norm.reverse()
        result.write("positions: " + ' '.join(str(e) for e in l_posn_norm) + '\n')
        for indv_num in range(1,len(d_hap_exon)+1):
            result.write(d_hap_exon[indv_num][::-1] + '\n')
    return ("written to file")

for exon in l_exons:
    print (exon)
    #exon part
    f_vcf = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/" + subfolder + "/" + s_yri_folder + "/VCF_FILTERED/" + exon + "_exon.vcf", 'r')
    t_exon = convert_vcf_to_ms(f_vcf)
    f_vcf.close()
    d_hap_exon = t_exon[0]
    l_posn_norm = normalize_positions(t_exon[1], int(d_exon_start[exon]), int(d_exon_end[exon]), "forward")
    result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/" + subfolder + "/" + s_yri_folder + "/MS_FILTERED/" + exon + "_exon.ms", 'w+')
    write_ms_file(result, l_posn_norm, d_hap_exon, "forward")
    result.close()
    
    #5prime intergenic:
    f_vcf = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/" + subfolder + "/" + s_yri_folder + "/VCF_FILTERED/" + exon + "_5p.vcf", 'r')
    t_5p = convert_vcf_to_ms(f_vcf)
    f_vcf.close()
    d_hap_5p = t_5p[0]
    if d_strand[exon] == "+":
        s_direction = "forward"
        s_5p_end = int(d_exon_start[exon]) - 1
        s_5p_start = int(d_exon_start[exon]) - 4*int(d_numbp50[exon])
    elif d_strand[exon] == "-":
        s_direction = "reverse"
        s_5p_start = int(d_exon_end[exon]) + 1
        s_5p_end = int(d_exon_end[exon]) + 4*int(d_numbp50[exon])
    l_posn_norm = normalize_positions(t_5p[1], int(s_5p_start), int(s_5p_end), s_direction)
    result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/" + subfolder + "/" + s_yri_folder + "/MS_FILTERED/" + exon + "_5p.ms", 'w+')
    write_ms_file(result, l_posn_norm, d_hap_5p, s_direction)
    result.close()
    
    #3prime intergenic:
    f_vcf = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/" + subfolder + "/" + s_yri_folder + "/VCF_FILTERED/" + exon + "_3p.vcf", 'r')
    t_3p = convert_vcf_to_ms(f_vcf)
    f_vcf.close()
    d_hap_3p = t_3p[0]
    if d_strand[exon] == "+":
        s_direction = "reverse"
        s_3p_start = int(d_exon_end[exon]) + 1
        s_3p_end = int(d_exon_end[exon]) + 4*int(d_numbp50[exon])
    elif d_strand[exon] == "-":
        s_direction = "forward"
        s_3p_end = int(d_exon_start[exon]) - 1
        s_3p_start = int(d_exon_start[exon]) - 4*int(d_numbp50[exon])
    l_posn_norm = normalize_positions(t_3p[1], int(s_3p_start), int(s_3p_end), s_direction)
    result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/" + subfolder + "/" + s_yri_folder + "/MS_FILTERED/" + exon + "_3p.ms", 'w+')
    write_ms_file(result, l_posn_norm, d_hap_3p, s_direction)
    result.close()


print("done")

