#This is to convert the vcf files into fasta files:
#Need to align your ms such that it's intergenic first and then the exon:
#The purpose is to be able to substract polymorphic sites from substituted sites.
#The result will have an "X" to indicate a SNP, "N" to indicate masked site (inaccessible), and "-" to indicate monomorphic site.

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

def get_positions_from_vcf(f_VCF):
    l_posns = []
    for line in f_VCF:
        line1 = line.strip('\n')
        if line1[0]!= "#":
            if "VT=SNP" in line1:#represents a SNP
                line2 = line1.split('\t')
                #check if the site is monomorphic:
                s_gts = ''.join(line2[9:])
                if "1" in s_gts:#site is polymorphic in YRI
                    l_posns.append(int(line2[1]))
    return(l_posns)

def get_filtered_sites(f_FILT):
    l_FILT = []
    for line in f_FILT:
        line1 = line.strip('\n')
        l_FILT.append(int(line1))
    return(l_FILT)

def get_sequence(l_FILT, l_POLY, posn_start, posn_end):
    SEQ = ""
    s_posn = int(posn_start)
    while s_posn <= int(posn_end):
        if s_posn in l_FILT:
            if s_posn in l_POLY:
                SEQ = SEQ + "X" #to indicate a SNP
            else:
                SEQ = SEQ + "-"
        else:
            SEQ = SEQ + "N" #masked
        s_posn = s_posn + 1
    return(SEQ)

def write_fasta_file(result, SEQ, direction):
    if direction == "forward":
        result.write(">" + '\n' + SEQ + '\n')
    elif direction == "reverse":
        result.write(">" + '\n' + SEQ[::-1] + '\n')
    return ("written to file")

for exon in l_exons:
    print (exon)
    #exon part
    #read filtered sites:
    f_filt = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/FILTERED_SITES/" + exon + "_exon.sites", 'r')
    l_filtered_exon = get_filtered_sites(f_filt)
    f_filt.close()
    if "exon1_" in exon:
        print(l_filtered_exon)
    #get polymorphic sites:
    f_vcf = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/" + subfolder + "/" + s_yri_folder + "/VCF_FILTERED/" + exon + "_exon.vcf", 'r')
    l_poly_exon = get_positions_from_vcf(f_vcf)
    f_vcf.close()
    #get final sequence
    s_seq = get_sequence(l_filtered_exon, l_poly_exon, int(d_exon_start[exon]), int(d_exon_end[exon]))
    result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/" + subfolder + "/" + s_yri_folder + "/FASTA_FILTERED/" + exon + "_exon.fa", 'w+')
    write_fasta_file(result, s_seq, "forward")
    result.close()
    
    #5prime intergenic:
    #read filtered sites:
    f_filt = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/FILTERED_SITES/" + exon + "_5p.sites", 'r')
    l_filtered_5p = get_filtered_sites(f_filt)
    f_filt.close()
    f_vcf = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/" + subfolder + "/" + s_yri_folder + "/VCF_FILTERED/" + exon + "_5p.vcf", 'r')
    l_poly_5p = get_positions_from_vcf(f_vcf)
    f_vcf.close()
    #write it
    if d_strand[exon] == "+":
        s_direction = "forward"
        s_5p_end = int(d_exon_start[exon]) - 1
        s_5p_start = int(d_exon_start[exon]) - 4*int(d_numbp50[exon])
    elif d_strand[exon] == "-":
        s_direction = "reverse"
        s_5p_start = int(d_exon_end[exon]) + 1
        s_5p_end = int(d_exon_end[exon]) + 4*int(d_numbp50[exon])
    s_seq = get_sequence(l_filtered_5p, l_poly_5p, int(s_5p_start), int(s_5p_end))
    result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/" + subfolder + "/" + s_yri_folder + "/FASTA_FILTERED/" + exon + "_5p.fa", 'w+')
    write_fasta_file(result, s_seq, s_direction)
    result.close()
    
    #3prime intergenic:
    #read filtered sites:
    f_filt = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/FILTERED_SITES/" + exon + "_3p.sites", 'r')
    l_filtered_3p = get_filtered_sites(f_filt)
    f_filt.close()
    f_vcf = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/" + subfolder + "/" + s_yri_folder + "/VCF_FILTERED/" + exon + "_3p.vcf", 'r')
    l_poly_3p = get_positions_from_vcf(f_vcf)
    f_vcf.close()
    #write it
    if d_strand[exon] == "+":
        s_direction = "reverse"
        s_3p_start = int(d_exon_end[exon]) + 1
        s_3p_end = int(d_exon_end[exon]) + 4*int(d_numbp50[exon])
    elif d_strand[exon] == "-":
        s_direction = "forward"
        s_3p_end = int(d_exon_start[exon]) - 1
        s_3p_start = int(d_exon_start[exon]) - 4*int(d_numbp50[exon])
    s_seq = get_sequence(l_filtered_3p, l_poly_3p, int(s_3p_start), int(s_3p_end))
    result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/" + subfolder + "/" + s_yri_folder + "/FASTA_FILTERED/" + exon + "_3p.fa", 'w+')
    write_fasta_file(result, s_seq, s_direction)
    result.close()


print("done")

