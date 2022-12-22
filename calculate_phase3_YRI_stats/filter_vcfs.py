#This is to filter vcf files according to the accessibility regions and phastcons elements:

import sys
s_yri_folder = sys.argv[1] #YRI/YRI_50

l_regions = ["exon", "5p", "3p"]

#get all exons:
f_list = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/FILTERED_SITES/num_of_filtered_sites.txt", 'r')
l_exons = []
for line in f_list:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] != "exon":
        l_exons.append(line2[0])
f_list.close()
#print l_exons

#get_all sites
for exon in l_exons:
    print (exon)
    for region in l_regions:
        d_sites = {}
        f_sites = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/FILTERED_SITES/" + exon + "_" + region + ".sites", 'r')
        for line in f_sites:
            line1 = line.strip('\n')
            d_sites[line1] = 1
        f_sites.close()
        #filter vcf file:
        f_vcf = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/VCF/" + exon + "_" + region + ".vcf", 'r')
        result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/VCF_FILTERED/" + exon + "_" + region + ".vcf", 'w+')
        for line in f_vcf:
            if line[0] == "#":
                result.write(line)
            else:
                line1 = line.strip('\n')
                line2 = line1.split('\t')
                if d_sites.get(line2[1], "NA") == 1:
                    result.write(line)
        f_vcf.close()
        result.close()

print ("done")



