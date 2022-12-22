#Basic stats, sliding window, ms format phase 3 data for YRI individuals
#only polymorphism based stats (no divergence)
#python statistics_bigwindow_pylibseq_SingExon_human_phase3.py -region exon

from __future__ import print_function
import libsequence
import sys
import pandas
import math
import argparse


#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-region', dest = 'region', action='store', nargs = 1, type = str, help = 'exon/5p/3p')
parser.add_argument('-yriFolder', dest = 'yriFolder', action='store', nargs = 1, type = str, help = 'YRI/YRI_50')
args = parser.parse_args()
s_region = args.region[0]
s_yri_folder = args.yriFolder[0]

#defined constants:

#get length of coding from file:
f_sing = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon2-6kb_CNE500_rec_numbp50_final_div.tab", 'r')
d_codingsize = {}
d_numbp50 = {}
d_filtered_sites = {}
l_exons = []
exonID = 1
for line in f_sing:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] == "exon":
        d_col = {}
        colnum = 0
        for x in line2:
            d_col[x] = colnum
            colnum = colnum + 1
    else:
        s_gene = line2[d_col["gene"]]
        s_exon = "exon" + str(exonID) + "_" + s_gene
        l_exons.append(s_exon)
        d_filtered_sites[s_exon] = {}
        d_codingsize[s_exon] = line2[d_col["exon_len"]]
        d_numbp50[s_exon] = line2[d_col["numbp50"]]
        exonID += 1
f_sing.close()

#get total number of usable sites from file:
f_num_sites = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/FILTERED_SITES/num_of_filtered_sites.txt", 'r')
for line in f_num_sites:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] != "exon":
        if line2[1] == s_region:
            exon = line2[0]
            window = line2[2]
            num_sites = line2[3]
            d_filtered_sites[exon][window] = num_sites
f_num_sites.close()

#read ms file:
def read_subset_ms(f_ms, start, end):
	l_Pos = [] #list of positions of SNPs
	#l_Genos = [] #list of alleles
	d_tmp = {}
	l_int_posn = []
	for line in f_ms:
		line1 = line.strip('\n')
		if "positions" in line1:
			line2 = line1.split()
			i = 0
			for x in line2:
				if "position" not in x:
					if (float(x) >= float(start)) and (float(x) <= float(end)):
						l_Pos.append(float(x))
						d_tmp[str(i)] = ""
						l_int_posn.append(i)
					i = i + 1
		elif "//" not in line and "segsites" not in line:
			for j in l_int_posn:
				d_tmp[str(j)] = d_tmp[str(j)] + line1[j]
	return (l_Pos, d_tmp, l_int_posn)



def avg_divergence_win(d_subs, start, end):
	s_sum = 0
	for posn in d_subs.keys():
		if float(posn) < float(end) and float(posn) > float(start):
			s_sum = s_sum + 1
	return s_sum

def print_5_dig(s_number):
    #print (s_number)
    if s_number == "NA" or s_number == "nan" or s_number == "inf":
        to_print = s_number
    else:
        #print(s_num)
        to_print = str("{:.5f}".format(s_number))
    return to_print

#name windows in which to calculate stats:
#four windows for intergenic regions: 4 window sof length numbp50
#single window for exonic regions
if s_region == "5p" or s_region == "3p":
    l_windows = ["linked4", "linked3", "linked2", "linked1"]
elif s_region == "exon":
    l_windows = ["func"]
result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/STATS_FILTERED/human_phase3_bigwindow_" + s_region + ".stats", 'w+')
result.write("exon" + '\t' + "WinType" + '\t' + "WinSize" + '\t' + "WinSize_filtered" + '\t' + "S" + '\t' + "thetapi" + '\t' + "thetaw" + '\t' + "thetah" + '\t' + "hprime" + '\t' + "tajimasd" +  '\t' + "numSing" + '\t' + "hapdiv" + '\t' + "rsq" + '\t' + "D" + '\t' + "Dprime" + '\n')


#go through all exons and read data into pylibseq format
s_absent = 0
for s_exon in l_exons:
    #define windows specific to the gene
    s_coding_size = int(d_codingsize[s_exon])
    s_numbp50 = int(d_numbp50[s_exon])
    d_start, d_end = {}, {}
    if s_region == "5p" or s_region == "3p":
        chr_len = 4.0*float(s_numbp50)
        d_start["linked4"], d_end["linked4"] = 0, float(s_numbp50-1)/float(chr_len-1)
        d_start["linked3"], d_end["linked3"] = float(s_numbp50)/float(chr_len), float((2*s_numbp50)-1)/float(chr_len-1)
        d_start["linked2"], d_end["linked2"] = float(2*s_numbp50)/float(chr_len), float((3*s_numbp50)-1)/float(chr_len-1)
        d_start["linked1"], d_end["linked1"] = float(3*s_numbp50)/float(chr_len), float((4*s_numbp50)-1)/float(chr_len-1)
    elif s_region == "exon":    
        chr_len = int(s_coding_size)
        d_start["func"], d_end["func"] = 0.0, 1.0
    
    print ("Reading files in: exon" + str(s_exon) + "for region: " + s_region)
    print ("bin starts and ends are:")
    print (d_start)
    print (d_end)
            
    #reading ms file:
    f_ms = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/" + s_yri_folder + "/MS_FILTERED/"+ s_exon + "_" + s_region + ".ms", 'r')
    d_Pos, d_tmp, d_intPos = {}, {}, {}
    for win in l_windows:
        t_ms = read_subset_ms(f_ms, d_start[win], d_end[win])
        d_Pos[win] = t_ms[0]
        d_tmp[win] = t_ms[1]
        d_intPos[win] = t_ms[2]
        f_ms.seek(0)
    f_ms.close()

    #Get sd -> another structure of pylibseq:
    d_sd = {}
    for win in l_windows:
        l_data = []
        l_genotypes = []
        j = 0
        for s_posn in d_intPos[win]:
            l_genotypes.append(d_tmp[win][str(s_posn)])
            t_tmp = (d_Pos[win][j], d_tmp[win][str(s_posn)])
            l_data.append(t_tmp)
            j = j + 1
        d_sd[win] = libsequence.SimData(l_data)
            
    #calculating stats:
    d_ps = {}
    for win in l_windows:
        d_ps[win] = libsequence.PolySIM(d_sd[win])
            
    #calculate LD stats:
    meanrsq, meanD, meanDprime = {}, {}, {}
    for win in l_windows:
        LD = libsequence.ld(d_sd[win])
        LDstats = pandas.DataFrame(LD)
        if len(LDstats) > 0:
            meanrsq[win] = sum(LDstats['rsq'])/len(LDstats['rsq'])
            meanD[win] = sum(LDstats['D'])/len(LDstats['D'])
            meanDprime[win] = sum(LDstats['Dprime'])/len(LDstats['Dprime'])
        else:
            meanrsq[win], meanD[win], meanDprime[win] = "NA", "NA", "NA"
            
    #write out stats:
    for win in l_windows:
        result.write(s_exon + '\t' + win + '\t' + str(round((d_end[win]-d_start[win])*float(chr_len))+1) + '\t' + d_filtered_sites[s_exon][win] + '\t' + str(d_ps[win].numpoly()) + '\t' + print_5_dig(d_ps[win].thetapi()) + '\t' + print_5_dig(d_ps[win].thetaw()) + '\t' + print_5_dig(d_ps[win].thetah()) + '\t' + print_5_dig(d_ps[win].hprime()) + '\t' + print_5_dig(d_ps[win].tajimasd()) + '\t' + str(d_ps[win].numexternalmutations()) + '\t' + print_5_dig(d_ps[win].hapdiv()))
        #write out LD stats:
        result.write('\t' + print_5_dig(meanrsq[win]) + '\t' + print_5_dig(meanD[win]) + '\t' + print_5_dig(meanDprime[win]) + '\n')
        #write out divergence:
        #result.write(print_5_dig(avg_divergence_win(d_subs[win], d_start[win], d_end[win])) + '\n')

result.close()
print ("Finished")






