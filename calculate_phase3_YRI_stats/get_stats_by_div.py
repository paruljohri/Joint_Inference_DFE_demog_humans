#This is to sort the exons by div and rec:

import sys

region = sys.argv[1] #5p/exon

l_categories = ["low_div", "high_div"]

d_category = {}

for cat in l_categories:
    f_exons = open("/home/pjohri1/BgsDfeDemo_Human/Humans/single_exon_" + cat + ".txt", 'r')
    for line in f_exons:
        line1 = line.strip('\n')
        line2 = line1.split('\t')
        if line2[0] != "exon":
            d_category[line2[0]] = cat
    f_exons.close()

#get divergence:
#5p:
f_div = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/divergence/YRI_50/divergence_" + region + "_filtered.txt", 'r')

result1 = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/divergence/YRI_50/by_category/divergence_" + region + "_filtered_low_div.txt", 'w+')
result2 = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/divergence/YRI_50/by_category/divergence_" + region + "_filtered_high_div.txt", 'w+')
d_files = {}
d_files["low_div"] = result1
d_files["high_div"] = result2

for line in f_div:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] == "exon":
        for cat in l_categories:
            d_files[cat].write(line)
    else:
        s_exon = line2[0].replace("_" + region, "")
        if d_category.get(s_exon, "NA") != "NA":
            d_files[d_category[s_exon]].write(line)
f_div.close()

for cat in l_categories:
    d_files[cat].close()

#write out polymorphism stats:
f_stats = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/YRI_50/STATS_FILTERED/human_phase3_bigwindow_" + region + ".stats", 'r')

result1 = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/YRI_50/STATS_FILTERED/by_category/human_phase3_bigwindow_" + region + "_low_div.stats", 'w+')
result2 = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/YRI_50/STATS_FILTERED/by_category/human_phase3_bigwindow_" + region + "_high_div.stats", 'w+')

d_files = {}                                                                    
d_files["low_div"] = result1                                            
d_files["high_div"] = result2

for line in f_stats:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] == "exon":
        for cat in l_categories:
            d_files[cat].write(line)
    else:
        s_exon = line2[0]
        if d_category.get(s_exon, "NA") != "NA":
            d_files[d_category[s_exon]].write(line)
f_stats.close()

for cat in l_categories:
    d_files[cat].close()

print("done")
