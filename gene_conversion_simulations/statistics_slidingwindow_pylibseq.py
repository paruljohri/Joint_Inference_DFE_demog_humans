#Basic stats, sliding window, SLIm ms output
#adding divergence to this
#python statistics_slidingwindow_pylibseq.py -winSize 500 -stepSize 500 -folder geneconv_humans -simID 3 -noncodingLen 50000 -codingLen 4000
from __future__ import print_function
#from libsequence.polytable import SimData
#from libsequence.summstats import PolySIM
#from libsequence.windows import Windows
#from libsequence.summstats import ld
import libsequence
import sys
import pandas
import math
import argparse


#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-winSize', dest = 'winSize', action='store', nargs = 1, default = 100, type = int, help = 'size of each sliding window in bp')#500 bp for small, 5000bp for big
parser.add_argument('-stepSize', dest = 'stepSize', action='store', nargs = 1, default = 100, type = int, help = 'size of step size in bp')#250 bp for small, 5000 bp for big
parser.add_argument('-folder', dest = 'FolderName', action='store', nargs = 1, type = str, help = 'the name of the folder or simulation to run')
parser.add_argument('-simID', dest = 'simulationID', action='store', nargs = 1, type = str, help = 'the name of the subfolder or simulation to run')
parser.add_argument('-noncodingLen', dest = 'noncodingLen', action='store', nargs = 1, type = int, help = 'noncoding length in bp')
parser.add_argument('-codingLen', dest = 'codingLen', action='store', nargs = 1, type = int, help = 'length in bp of coding region simulated')#Length of coding region simulated
args = parser.parse_args()
noncoding_size =  args.noncodingLen[0]
coding_size =  args.codingLen[0]
chr_len = noncoding_size + coding_size
win_size = args.winSize[0]/float(chr_len)
step_size = args.stepSize[0]/float(chr_len)
folder = args.FolderName[0]
simID = args.simulationID[0]

#DEFINE!
num_reps = 1000

def read_fixed_mutations(f_fixed):
    d_subs = {}
    for line in f_fixed:
        line1 = line.strip('\n')
        line2 = line1.split()
        if line1[0]!="#" and line2[0]!="Mutations:":
            posn = float(line2[3])/float(chr_len)
            if "grow10" in folder  or "red10" in folder:
                num_gen = line2.pop()
                if int(num_gen) >= 50000:
                    d_subs[posn] = d_subs.get(posn, 0) + 1
            else:
                d_subs[posn] = d_subs.get(posn, 0) + 1
    return d_subs #return a dictionary with base positions as keys and number of fixed substitutions as values

def avg_divergence_win(d_subs, start, end):
	s_sum = 0
	for posn in d_subs.keys():
		if float(posn) <= end and float(posn) > start:
			s_sum = s_sum + 1
			
	return s_sum

#result files:
result = open("/scratch/pjohri1/BgsDfeDemo_Human/" + folder + "/sim" + str(simID) + "_win" + str(args.winSize[0]) + "_step" + str(args.stepSize[0]) + ".stats", 'w+')
result.write("simnum" + '\t' + "posn" + '\t' + "thetapi" + '\t' + "thetaw" + '\t' + "thetah" + '\t' + "hprime" + '\t' + "tajimasd" +  '\t' + "numSing" + '\t' + "hapdiv" + '\t' + "rsq" + '\t' + "D" + '\t' + "Dprime" + '\t' + "div" + '\n')

#go through all simulation replicates and read data into pylibseq format
#addin the option of ignoring some files if they don't exist
repID = 1
s_absent = 0
while repID <= num_reps:
    print ("Reading file:" + str(repID))
    #try:
    if repID > 0:
        f_ms = open("/scratch/pjohri1/BgsDfeDemo_Human/" + folder + "/sim" + str(simID) + "/sim" + str(simID) + "_rep" + str(repID) + ".ms", 'r')
        f_subs = open("/scratch/pjohri1/BgsDfeDemo_Human/" + folder + "/sim" + str(simID) + "/sim" + str(simID) + "_rep" + str(repID) + ".fixed", 'r')
        d_subs = read_fixed_mutations(f_subs)
        l_Pos = [] #list of positions of SNPs
        l_Genos = [] #list of alleles
        d_tmp = {}
        for line in f_ms:
            line1 = line.strip('\n')
            if "positions" in line1:
                line2 = line1.split()
                i = 0
                for x in line2:
                    if "position" not in x:
                        l_Pos.append(float(x))
                        d_tmp[str(i)] = ""
                        i = i + 1
            elif "//" not in line and "segsites" not in line:
                #print (d_tmp)
                i = 0
                while i < len(line1):
                    d_tmp[str(i)] = d_tmp[str(i)] + line1[i]
                    i = i + 1
		#print (d_tmp)
        l_data = []
        i = 0
        while i < len(l_Pos):
            l_Genos.append(d_tmp[str(i)])
            t_tmp = (l_Pos[i], d_tmp[str(i)])
            l_data.append(t_tmp)
            i = i + 1
	#print (l_Pos)
	#print (l_Genos)


		#assign object
        sd = libsequence.SimData(l_data)
		#sd.assign(l_Pos[10:100],l_Genos[10:100])

		#define sliding windows:
        w = libsequence.Windows(sd,window_size=win_size,step_len=step_size,starting_pos=0.0,ending_pos=1.0)
		#chromosome length = 30kb, window size = 5 kb
        num_win = len(w)

		#calculate summary statistic in sliding window:
        print ("calculating stats in windows")
        win_name = 1
        for i in range(len(w)):
            wi = w[i]
            #print (wi)
            pswi = libsequence.PolySIM(wi)
            result.write("rep" + str(repID) + '\t' + str(win_name) + '\t' + str(pswi.thetapi()) + '\t' + str(pswi.thetaw()) + '\t' + str(pswi.thetah()) + '\t' + str(pswi.hprime()) + '\t' + str(pswi.tajimasd()) + '\t' + str(pswi.numexternalmutations()) + '\t' + str(pswi.hapdiv()) + '\t')
			#read data to calculate LD based stats:
            
            if len(wi.pos()) >= 5: #These are pairwise stats. If only 1 site exists, it'll show an error.
                #print (i)
                LD_tmp = libsequence.ld(wi)
                LDstats = pandas.DataFrame(LD_tmp)
                meanrsq = sum(LDstats['rsq'])/len(LDstats['rsq'])
                meanD = sum(LDstats['D'])/len(LDstats['D'])
                meanDprime = sum(LDstats['Dprime'])/len(LDstats['Dprime'])
                result.write(str(meanrsq) + '\t' + str(meanD) + '\t' + str(meanDprime) + '\t')
            else:
                result.write("NA" + '\t' + "NA" + '\t' + "NA" + '\t') 
        	#divergence:
            s_start = (i)*(1.0/float(num_win))
            s_end = s_start + win_size
            result.write(str(avg_divergence_win(d_subs, s_start, s_end)) + '\n')
            win_name = win_name + 1
    #except:
    else:
        s_absent = s_absent + 1
        print ("This file does not exist or cannot be read or is empty")
	
    repID = repID + 1

print ("Number of files not read:" + '\t' + str(s_absent))
print ("Finished")






