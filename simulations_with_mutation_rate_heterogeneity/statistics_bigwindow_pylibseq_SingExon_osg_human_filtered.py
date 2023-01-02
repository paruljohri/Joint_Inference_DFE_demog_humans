#Basic stats, sliding window, SLIm ms output
#adding divergence to this
#python statistics_bigwindow_pylibseq_SingExon_osg_human.py -folder demo_dfe_SingExon_human_v2 -simID 1 -numGenes 465 -numRep 1
from __future__ import print_function
import libsequence
import sys
import pandas
import math
import argparse
import os

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-folder', dest = 'FolderName', action='store', nargs = 1, type = str, help = 'the name of the folder or simulation to run')
parser.add_argument('-outFolder', dest = 'outFolder', action='store', nargs = 1, type = str, help = 'the name of the folder in which output is placed')
parser.add_argument('-interRegion', dest = 'interRegion', action='store', nargs = 1, type = str, help = '5p/3p; the intergenic region on which filtering is based')
parser.add_argument('-simID', dest = 'simulationID', action='store', nargs = 1, type = str, help = 'the name of the subfolder or simulation to run')
parser.add_argument('-numGenes', dest = 'numGenes', action='store', nargs = 1, type = int, help = 'number of genes')
parser.add_argument('-numRep', dest = 'numRep', action='store', nargs = 1, type = int, choices = range(1,10001), help = 'number of replicates of each gene')
args = parser.parse_args()
folder = args.FolderName[0]
out_folder = args.outFolder[0]
inter_region = args.interRegion[0]
simID = args.simulationID[0]
subfolder = "sim" + str(simID)
num_genes = args.numGenes[0]
num_rep = args.numRep[0]
#defined constants:

#get length of coding from file:
f_sing = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon2-6kb_CNE500_rec_numbp50_final_div.tab", 'r')
d_codingsize = {}
d_numbp50 = {}
geneID = 0
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
        geneID = geneID + 1
        d_codingsize["gene" + str(geneID)] = line2[d_col["exon_len"]]
        d_numbp50["gene" + str(geneID)] = line2[d_col["numbp50"]]
f_sing.close()


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

#get burn-in period from file:
f_par = open("/home/pjohri1/BgsDfeDemo_Human/programs_model_violations/model_violations_parameters.txt", 'r')
d_col = {}
for line in f_par:
    line1 = line.strip('\n')
    line2 = line1.split()
    if line2[0] == "simID":
        col = 0
        for x in line2:
            d_col[x] = col
            col = col + 1
    elif line2[0] == "sim" + str(simID):
        gen_burnin = 10*int(line2[d_col["Na_scaled"]])
f_par.close()
print ("gen_burnin: " + str(gen_burnin))

#read fixed mutations: output of Slim:
def read_subset_fixed_mutations(f_fixed, start, end):
    d_subs = {}
    for line in f_fixed:
        line1 = line.strip('\n')
        line2 = line1.split()
        if line1[0]!="#" and line2[0]!="Mutations:":
            posn = float(line2[3])/float(chr_len-1)
            num_gen = line2[8]
            #print (num_gen)
            if int(num_gen) >= gen_burnin:#include this mutation only if it fixed after burnin period- 10Na.
                #print (line2)
                if posn >= float(start) and posn <= float(end):
                    d_subs[posn] = d_subs.get(posn, 0) + 1
    return d_subs #return a dictionary with base positions as keys and number of fixed substitutions as values

def avg_divergence_win(d_subs, start, end):
	s_sum = 0
	for posn in d_subs.keys():
		if float(posn) <= float(end) and float(posn) >= float(start):
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

#get the total number of sites (after filtering):
os.system("cp /home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/FILTERED_SITES/num_of_filtered_sites.txt /scratch/pjohri1/BgsDfeDemo_Human/" + folder + "/" + subfolder)
f_filt = open("/scratch/pjohri1/BgsDfeDemo_Human/" + folder + "/" + subfolder + "/num_of_filtered_sites.txt", 'r')
d_filt_sites = {}                                                               
for line in f_filt:
    line1 = line.strip('\n')                                                    
    line2 = line1.split('\t')                                                   
    if line2[0] != "exon":
        if line2[1] == "exon" or line2[1] == inter_region:
            exon_num = line2[0].split('_')[0]
            try:                                                                    
                d_filt_sites[exon_num.replace("exon", "gene")][line2[2]] = line2[3] 
            except:                                                                 
                d_filt_sites[exon_num.replace("exon", "gene")] = {}                 
                d_filt_sites[exon_num.replace("exon", "gene")][line2[2]] = line2[3] 
f_filt.close()                                                                  
print(d_filt_sites)

#name windows in which to calculate stats:
#five windows: 4 window sof length numbp50 and then the exon
l_windows = ["linked4", "linked3", "linked2", "linked1", "func"]

result = open("/scratch/pjohri1/BgsDfeDemo_Human/" + out_folder + "/" + subfolder + "_bigwindow.stats", 'w+')
result.write("geneID" + '\t' + "WinType" + '\t' + "WinSize" + '\t' + "WinSize_filtered"+ '\t' + "thetapi" + '\t' + "thetaw" + '\t' + "thetah" + '\t' + "hprime" + '\t' + "tajimasd" +  '\t' + "numSing" + '\t' + "hapdiv" + '\t' + "rsq" + '\t' + "D" + '\t' + "Dprime" + '\t' + "div" + '\n')


#go through all simulation replicates and read data into pylibseq format
#addin the option of ignoring some files if they don't exist
geneID = 1
s_absent = 0
while geneID <= num_genes:
    #define windows specific to the gene
    s_coding_size = int(d_codingsize["gene" + str(geneID)])
    s_numbp50 = int(d_numbp50["gene" + str(geneID)])
    chr_len = 4.0*float(s_numbp50) + int(s_coding_size)
    d_start, d_end = {}, {}
    d_start["linked4"], d_end["linked4"] = 0, float(s_numbp50-1)/float(chr_len-1)
    d_start["linked3"], d_end["linked3"] = float(s_numbp50)/float(chr_len-1), float((2*s_numbp50)-1)/float(chr_len-1)
    d_start["linked2"], d_end["linked2"] = float(2*s_numbp50)/float(chr_len-1), float((3*s_numbp50)-1)/float(chr_len-1)
    d_start["linked1"], d_end["linked1"] = float(3*s_numbp50)/float(chr_len-1), float((4*s_numbp50)-1)/float(chr_len-1)
    d_start["func"], d_end["func"] = float(4*s_numbp50)/float(chr_len-1), 1.0
    
    print ("Reading files in: gene" + str(geneID))
    print ("bin starts and ends are:")
    print (d_start)
    print (d_end)
    
    repID = 1
    while repID <= num_rep:
        if repID > 0:
        #try:
            print ("Reading file: sim" + str(simID) + "_gene" + str(geneID) + "_rep" + str(repID))
            f_subs = open("/scratch/pjohri1/BgsDfeDemo_Human/" + folder + "/" + subfolder + "/sim" + str(simID) + "_gene" + str(geneID) + "_rep" + str(repID) + ".fixed", 'r')
            d_subs = {}
            for win in l_windows:
                d_subs[win] = read_subset_fixed_mutations(f_subs, d_start[win], d_end[win])
                f_subs.seek(0)
                #print (len(d_subs[win]))
            f_subs.close()
            
            #reading ms file:
            f_ms = open("/scratch/pjohri1/BgsDfeDemo_Human/" + folder + "/" + subfolder + "/sim" + str(simID) + "_gene" + str(geneID) + "_rep" + str(repID) + ".ms", 'r')
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
                result.write("gene" + str(geneID) + "_rep" + str(repID) + '\t' + win + '\t' + str(round((d_end[win]-d_start[win])*float(chr_len))) + '\t' + str(d_filt_sites["gene" + str(geneID)][win]) + '\t' + print_5_dig(d_ps[win].thetapi()) + '\t' + print_5_dig(d_ps[win].thetaw()) + '\t' + print_5_dig(d_ps[win].thetah()) + '\t' + print_5_dig(d_ps[win].hprime()) + '\t' + print_5_dig(d_ps[win].tajimasd()) + '\t' + str(d_ps[win].numexternalmutations()) + '\t' + print_5_dig(d_ps[win].hapdiv()))
                #write out LD stats:
                result.write('\t' + print_5_dig(meanrsq[win]) + '\t' + print_5_dig(meanD[win]) + '\t' + print_5_dig(meanDprime[win]) + '\t')
                #write out divergence:
                result.write(print_5_dig(avg_divergence_win(d_subs[win], d_start[win], d_end[win])) + '\n')
        else:
        #except:
            s_absent = s_absent + 1
            print ("This file does not exist or cannot be read or is empty")
        repID = repID + 1
    geneID = geneID + 1

result.close()
print ("Number of files not read:" + '\t' + str(s_absent))
print ("Finished")






