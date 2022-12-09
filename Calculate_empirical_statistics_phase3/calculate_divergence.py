#This is to calcualte divergence per site, by comparing hte refernece genome to the ancestral genome.
#We will try using both the hg19 genome and hte grch37 genome.

import sys
s_popn = sys.argv[1] #YRI/YRI_50
s_region = sys.argv[2] #exon/5p/3p

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

def calculate_divergence(SEQ1, SEQ2, SEQ_SNPS):
    s_diff, s_tot = 0, 0
    if len(SEQ1) == len(SEQ2) and len(SEQ1) == len(SEQ_SNPS):
        for i in range(0, len(SEQ1), 1):
            if SEQ1[i] not in [".", "-", "N", "n"] and SEQ2[i] not in [".", "-", "N", "n"]:#both sequences have a known base
                if SEQ_SNPS[i] != "X" and SEQ_SNPS[i] != "N":#the site is fixed
                    s_tot += 1
                    if SEQ1[i].upper() != SEQ2[i].upper():
                        s_diff += 1
    else:
        print("check, length of sequences or SNPs is not the same")
    if s_tot == 0:
        return("NA")
    else:
        return(float(s_diff)/float(s_tot))

#open file to write divergence values:
result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/divergence/" + s_popn + "/divergence_" + s_region + "_filtered.txt", 'w+')
if s_region == "5p" or s_region == "3p":
    result.write("exon" + '\t' + "linked4_hg19-anc" + '\t' + "linked3_hg19-anc" + '\t' + "linked2_hg19-anc" + '\t' + "linked1_hg19-anc" + '\t' + "linked4_grch37-anc" + '\t' + "linked3_grch37-anc" + '\t' + "linked2_grch37-anc" + '\t' + "linked1_grch37-anc" + '\n')
elif s_region == "exon":
    result.write("exon" + '\t' + "func_hg19-anc" + '\t' + "func_grch37-anc" + '\n')

#read alignement files and store divergence values:
d_div_hg19, d_div_grch37 = {}, {}
chr_num = 1
while chr_num <= 22: #should be 22
    f_anc = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/divergence/ANC/" + s_region + "/chr" + str(chr_num) + ".fa", 'r')
    f_hg19 = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/divergence/HG19/" + s_region + "/chr" + str(chr_num) + ".fa", 'r')
    f_grch37 = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/divergence/GRCH37/" + s_region + "/chr" + str(chr_num) + ".fa", 'r')
    
    d_anc = read_fasta(f_anc)
    d_hg19 = read_fasta(f_hg19)
    d_grch37 = read_fasta(f_grch37)

    f_anc.close()
    f_hg19.close()
    f_grch37.close()
    
    #get divergence
    for x in d_anc.keys():
        s_exon = x.replace("_anc", "")
        f_poly = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_popn + "/FASTA_FILTERED/" + s_exon + ".fa", 'r')
        d_poly = read_fasta(f_poly)
        f_poly.close()
        #print(d_poly)
        if s_region == "exon":
            d_div_hg19[s_exon] = calculate_divergence(d_anc[x], d_hg19[s_exon + "_hg19"], d_poly[""])
            d_div_grch37[s_exon] = calculate_divergence(d_anc[x], d_grch37[s_exon + "_grch37"], d_poly[""])
            result.write(s_exon + '\t' + str(d_div_hg19[s_exon]) + '\t' + str(d_div_grch37[s_exon]) + '\n')
        elif s_region == "5p" or s_region == "3p":
            result.write(s_exon)
            win_size = round(len(d_anc[x])/4.0)
            print (str(len(d_anc[x])) + '\t' + str(win_size))
            i = 0
            while i < 4:
                s_div_hg19 = calculate_divergence(d_anc[x][i*win_size:(i+1)*win_size], d_hg19[s_exon + "_hg19"][i*win_size:(i+1)*win_size], d_poly[""][i*win_size:(i+1)*win_size])
                result.write('\t' + str(s_div_hg19))
                i += 1
            i = 0
            while i < 4:
                s_div_grch37 = calculate_divergence(d_anc[x][i*win_size:(i+1)*win_size], d_grch37[s_exon + "_grch37"][i*win_size:(i+1)*win_size], d_poly[""][i*win_size:(i+1)*win_size])
                result.write('\t' + str(s_div_grch37))
                i += 1
            result.write('\n')
    chr_num += 1

result.close()
print ("done")


