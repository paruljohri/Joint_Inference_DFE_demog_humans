#This is to select the list of exons that should finally qualify for analysis:
#Also add filters for average recombination rate: >= 0.5cM/Mb and <= 10 cM/Mb

import sys

inter_cutoff = 20000


result50 = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon2-6kb_CNE500_rec_numbp50_final.tab", 'w+')
result_fixed = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon2-6kb_CNE500_rec_inter20000_final.tab", 'w+')

d_col = {}

f_tab = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon1-6kb_CNE500_rec_numbp50.tab", 'r')
for line in f_tab:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] == "exon":
        result50.write(line1 + '\t' + "length_5p" + '\t' + "length_3p" + '\n')
        result_fixed.write(line1 + '\t' + "length_5p" + '\t' + "length_3p" + '\n')
        col = 0
        for x in line2:
            d_col[x] = col
            col += 1
    else:
        print (line2[0])
        if "begin" not in line and "end" not in line:
            if int(line2[d_col["exon_len"]]) >= 2000 and int(line2[d_col["exon_len"]]) <= 6000 and line2[d_col["avg_rec_rate"]] != "NA":#exons are of the size we want
                if float(line2[d_col["avg_rec_rate"]]) >= 0.5 and float(line2[d_col["avg_rec_rate"]]) <= 10.0:
                    if line2[d_col["strand"]] == "+":
                        if line2[d_col["5p_elements"]] == "none":
                            inter_5p = line2[d_col["intergenic_5p_len"]]
                        else:
                            inter_5p = abs(int(line2[d_col["intergenic_5p_end"]]) - int(line2[d_col["5p_end"]]))
                        if line2[d_col["3p_elements"]] == "none":
                            inter_3p = line2[d_col["intergenic_3p_len"]]
                        else:
                            inter_3p = abs(int(line2[d_col["intergenic_3p_start"]]) - int(line2[d_col["3p_end"]]))
                        if int(inter_5p) >= inter_cutoff and int(inter_3p) >= inter_cutoff:
                            result_fixed.write(line1 + '\t' + str(inter_5p) + '\t' + str(inter_3p) + '\n')
                        if line2[d_col["numbp50"]] != "NA":
                            if int(inter_5p)>= int(line2[d_col["numbp50"]])*4 and int(inter_3p) >= int(line2[d_col["numbp50"]])*4:
                                result50.write(line1 + '\t' + str(inter_5p) + '\t' + str(inter_3p) + '\n')
                    elif line2[d_col["strand"]] == "-":
                        if line2[d_col["5p_elements"]] == "none":
                            inter_5p = line2[d_col["intergenic_5p_len"]]
                        else:
                            inter_5p = abs(int(line2[d_col["intergenic_5p_start"]]) - int(line2[d_col["5p_end"]]))
                        if line2[d_col["3p_elements"]] == "none":
                            inter_3p = line2[d_col["intergenic_3p_len"]]
                        else:
                            inter_3p = abs(int(line2[d_col["intergenic_3p_end"]]) - int(line2[d_col["3p_end"]]))
                        if int(inter_5p) >= inter_cutoff and int(inter_3p) >= inter_cutoff:
                            result_fixed.write(line1 + '\t' + str(inter_5p) + '\t' + str(inter_3p) + '\n')
                        if line2[d_col["numbp50"]] != "NA":
                            if int(inter_5p) >= int(line2[d_col["numbp50"]])*4 and int(inter_3p) >= int(line2[d_col["numbp50"]])*4:
                                result50.write(line1 + '\t' + str(inter_5p) + '\t' + str(inter_3p) + '\n')
f_tab.close()
result50.close()
result_fixed.close()
print ("Finished")

