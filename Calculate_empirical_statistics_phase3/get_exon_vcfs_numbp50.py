#This is to make .ms files from vcf files for the selected exons:
#This program does the same thing as _v1 but is much more efficient.
#In this script, we get stats from 4*numbp50 from both sides of the intergenic.

import sys
import os

s_yri_folder = sys.argv[1] #YRI/YRI_50

s_chrom = "chr" + str(sys.argv[2]) #chromosome number

#read in the individuals from YRI
d_yri = {}
if s_yri_folder == "YRI_50":
    f_yri = open("/scratch/pjohri1/BgsDfeDemo_Human/HumanData/Variants/YRI_individuals_50.txt",'r')
else:
    f_yri = open("/scratch/pjohri1/BgsDfeDemo_Human/HumanData/Variants/YRI_individuals.txt",'r')
for line in f_yri:
    line1 = line.strip('\n')
    d_yri[line1] = 1
f_yri.close()

#get chromosome names
d_chrom = {}
f_info = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_info.txt", 'r')
for line in f_info:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    d_chrom[line2[6]] = line2[9]
f_info.close()
print (d_chrom)

#get all positions you want
f_exons = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon2-6kb_CNE500_rec_numbp50_final_div.tab", 'r')
d_col = {}
l_exons = []
exon_num = 1
d_site_type = {}
for line in f_exons:
    lineA = line.strip('\n')
    lineB = lineA.split('\t')
    if lineB[0] == "exon":
        col_num = 0
        for x in lineB:
            d_col[x] = col_num
            col_num += 1
    else:
        s_chr = d_chrom[lineB[d_col["chrom"]]]
        if s_chr == s_chrom:
            s_exon_start = int(lineB[d_col["exon_start"]])
            s_exon_end = int(lineB[d_col["exon_end"]])
            s_numbp50 = int(lineB[d_col["numbp50"]])
            s_strand = lineB[d_col["strand"]]
            s_gene = lineB[d_col["gene"]]
            s_exon = "exon" + str(exon_num) + "_" + s_gene
            print("gene: " + s_exon)
            l_exons.append(s_exon)

            #get positions of exons vs intergenic:
            s_posn = s_exon_start
            while s_posn <= s_exon_end:
                d_site_type[s_chr + "_" + str(s_posn)] = s_exon + "_exon"
                s_posn += 1
            #going towards the left of the exonic sequence:
            s_posn = s_exon_start - 1
            inter_len = 0
            while inter_len <= 4*s_numbp50:
                if s_strand == "+":
                    d_site_type[s_chr + "_" + str(s_posn)] = s_exon + "_5p"
                elif s_strand == "-":
                    d_site_type[s_chr + "_" + str(s_posn)] = s_exon + "_3p"
                s_posn = s_posn - 1
                inter_len += 1
            #going towards the right of the exonic sequence:
            s_posn = s_exon_end + 1
            inter_len = 0
            while inter_len <= 4*s_numbp50:
                if s_strand == "+":
                    d_site_type[s_chr + "_" + str(s_posn)] = s_exon + "_3p"
                elif s_strand == "-":
                    d_site_type[s_chr + "_" + str(s_posn)] = s_exon + "_5p"
                s_posn += 1
                inter_len += 1
        
            #Create new files:
            f_exon = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/VCF/" + s_exon + "_exon.vcf", 'w+')
            f_5p = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/VCF/" + s_exon + "_5p.vcf", 'w+')
            f_3p = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/VCF/" + s_exon + "_3p.vcf", 'w+')
            f_exon.close()
            f_5p.close()
            f_3p.close()
        exon_num += 1
f_exons.close()
#print(d_site_type)


#uncompress the variant vcf file:
print("uncompressing vcf file chr" + s_chrom)
os.system("gunzip /scratch/pjohri1/BgsDfeDemo_Human/HumanData/Variants/ALL." + s_chrom + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")

#write sites into files:
print("going through the vcf file")
l_yri_posns = []
f_vcf = open("/scratch/pjohri1/BgsDfeDemo_Human/HumanData/Variants/ALL." + s_chrom + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf", 'r')
for line in f_vcf:
    line1 = line.strip('\n')
    if line1[0:2] != "##":
        line2 = line1.split('\t')
        if line2[0] == "#CHROM":
            #read in the column numbers of relevant individuals:
            col_num = 9
            for x in line2[9:]:
                if d_yri.get(x, "NA") == 1:
                    l_yri_posns.append(col_num)
                col_num += 1
            
            #write the column names into vcf files:
            print("creating new files and writing column names")
            for gene in l_exons:
                for region in ["exon", "5p", "3p"]:
                    result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/VCF/" + gene + "_" + region + ".vcf", 'a')
                    result.write('\t'.join(line2[0:9]))
                    for posn in l_yri_posns:
                        result.write('\t' + line2[posn])
                    result.write('\n')
                    result.close()
            print ("done writing column names")
        else:
            #identify which regions each site belongs to and write it in the relevant file:
            s_site = s_chrom + "_" + line2[1]
            s_site_type = d_site_type.get(s_site, "NA")
            #print(s_site + '\t' + s_site_type)
            if s_site_type != "NA":
                #print(s_site_type)
                #print(line)
                result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/VCF/" + s_site_type + ".vcf", 'a')
                result.write('\t'.join(line2[0:9]))
                for posn in l_yri_posns:
                    result.write('\t' + line2[posn])
                result.write('\n')
                result.close()
f_vcf.close()    

#compress the variant vcf file:  
print ("compressing vcf file again")                                     
os.system("gzip -v /scratch/pjohri1/BgsDfeDemo_Human/HumanData/Variants/ALL." + s_chrom + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf")


print("done")


