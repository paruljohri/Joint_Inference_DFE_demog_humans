#This is to get the final set of polymorphism and divergence stats in abc format:
#These currently include only hg19 divergence values. Can be changed.
#These currently include only windows linked1 and linked2. Can be changed.

import sys
s_yri_folder = sys.argv[1] #YRI/YRI_50
d_stats = {}

def read_bigwinsummary_file(FILE):
    l_windows = []
    for line in FILE:
        line1 = line.strip('\n')
        line2 = line1.split('\t')
        if "func" in line1 or "linked" in line1:
            for x in line2:
                l_windows.append(x)
        else:
            col = 0
            stat = line2[0]
            for x in line2[1:]:
                d_stats[l_windows[col] + "_" + stat] = x
                col += 1                
    return(d_stats)

def read_divergence_file(FILE):
    l_columns = []
    l_col_names = []
    for line in FILE:
        line1 = line.strip('\n')
        line2 = line1.split('\t')
        if "hg19" in line:
            col = 0
            for x in line2:
                x1 = x.replace("hg19.anc", "div")
                x2 = x1.replace("mean", "m")
                l_col_names.append(x2)
                if "hg19" in x:
                    l_columns.append(col)
                col += 1
        else:
            for col in l_columns:
                d_stats[l_col_names[col]] = line2[col]
    return(d_stats)

#read in all files:
f_stats_exon = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/STATS_FILTERED/human_phase3_bigwindow_exon.bigwinsummary", 'r')
f_stats_5p = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/STATS_FILTERED/human_phase3_bigwindow_5p.bigwinsummary", 'r')
f_stats_3p = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/STATS_FILTERED/human_phase3_bigwindow_3p.bigwinsummary", 'r')
f_div_exon = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/STATS_FILTERED/human_phase3_bigwindow_exon.divergence", 'r')
f_div_5p = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/STATS_FILTERED/human_phase3_bigwindow_5p.divergence", 'r')
f_div_3p = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/STATS_FILTERED/human_phase3_bigwindow_3p.divergence", 'r')


result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/STATS_FILTERED/human_phase3_abc.txt", 'w+')
                                                                                
l_stats = ["thetapi_m","thetaw_m","thetah_m", "hprime_m", "tajimasd_m", "numSing_m", "hapdiv_m", "rsq_m", "D_m", "Dprime_m", "div_m", "thetapi_sd", "thetaw_sd", "thetah_sd", "hprime_sd", "tajimasd_sd", "numSing_sd", "hapdiv_sd", "rsq_sd", "D_sd", "Dprime_sd", "div_sd"]
result.write("data")
for stat in l_stats:                                                            
    result.write('\t' + "func_" + stat + '\t' + "linked1_" + stat + '\t' + "linked2_" + stat)
result.write('\n') 

#5 prime:
d_stats = {}
read_bigwinsummary_file(f_stats_exon)
read_bigwinsummary_file(f_stats_5p)
read_divergence_file(f_div_exon)
read_divergence_file(f_div_5p)
f_stats_exon.seek(0)
f_stats_5p.close()
f_div_exon.seek(0)
f_div_5p.close()
#print(d_stats)
result.write("5prime")
for stat in l_stats:                                                            
    result.write('\t' + d_stats["func_" + stat] + '\t' + d_stats["linked1_" + stat] + '\t' + d_stats["linked2_" + stat])
result.write('\n')

#3prime:
d_stats = {}
read_bigwinsummary_file(f_stats_exon)
read_bigwinsummary_file(f_stats_3p)
read_divergence_file(f_div_exon)
read_divergence_file(f_div_3p)
f_stats_exon.close()
f_stats_3p.close()
f_div_exon.close()
f_div_3p.close()
result.write("3prime")
for stat in l_stats:
    result.write('\t' + d_stats["func_" + stat] + '\t' + d_stats["linked1_" + stat] + '\t' + d_stats["linked2_" + stat])
result.write('\n')

print ("done")




