#This is to select the list of exons that should qualify for intital screening:
#Length restriction on exons.
#Check if they are preotein-coding.

import sys

cutoff = 500

result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon1-6kb_CNE500.tab", 'w+')
chr_num = 1
d_col = {}
while chr_num <= 22:
    f_tab = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exons_firstelement_chr" + str(chr_num) + "_" + str(cutoff) + "bp.tab", 'r')
    for line in f_tab:
        line1 = line.strip('\n')
        line2 = line1.split('\t')
        if line2[0] == "exon":
            if chr_num == 1:
                result.write(line)
            col = 0
            for x in line2:
                d_col[x] = col
                col += 1
        else:
            print (line2[0])
            if "begin" not in line and "end" not in line:
                if line2[d_col["exon_type"]] == "protein_coding":
                    if int(line2[d_col["exon_len"]]) >= 1000 and int(line2[d_col["exon_len"]]) <= 6000:#exons are of the size we want
                        result.write(line)
    f_tab.close()
    chr_num += 1
result.close()

print ("Finished")

