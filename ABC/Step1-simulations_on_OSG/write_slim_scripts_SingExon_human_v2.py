#This is to perform 10 replicates of each exon:
#This program will create a temporary slim script and save it's output, and then delete the script.
#python write_slim_scripts_SingExon_human_v2.py -numRep 10 -simID 1
#go over starting ID to end ID. Check if a set has already been simulated. Then submit it if it hasn't been done.
#For reduced priors of 5-50k

from random import *
import os
import sys
import argparse
import math
parser = argparse.ArgumentParser(description='Information about number of replicates and cores used')
parser.add_argument('-numRep', dest = 'numReplicates', action='store', nargs = 1, default = 1, type = int, choices = range(1,10001), help = 'number of replicates of the same simulation')
parser.add_argument('-simID', dest = 'simID', action='store', nargs = 1, default = 1, type = int, help = 'the ID of the simulation to start with')

args = parser.parse_args()
num_rep = args.numReplicates[0]
simID = args.simID[0]

#Other constants:
num_gen_change = 200
scaling_factor = 1

#make a directory for results of this simualtion:
os.system("mkdir /home/pjohri/DEMO_DFE_SINGEXON_HUMAN_V2/sim" + str(simID))

#get parameters and ID of that simulation:
f_para = open("/home/pjohri/demo_dfe_SingExon_human_logunif_parameters_v2.txt", 'r')
d_col = {}
for line in f_para:
	line1 = line.strip('\n')
	line2 = line1.split('\t')
	if line2[0] == "simID":
		col = 0
		for x in line2:
			d_col[x] = col
			col = col + 1
	elif line2[0] == "sim" + str(simID):
		print (line)
		f0 = line2[d_col["f0"]]
		f1 = line2[d_col["f1"]]
		f2 = line2[d_col["f2"]]
		f3 = line2[d_col["f3"]]
		Na = int(line2[d_col["Na_scaled"]])
		Nc = int(line2[d_col["Nc_scaled"]])
		scaling_factor = line2[d_col["scaling_factor"]]
		g_factor = float(line2[d_col["g_factor"]])
f_para.close()

#re-write the main script and write out specific parameters in it with the simulation ID in it
f_exons = open("/home/pjohri/GRCh37_exon2-6kb_CNE500_rec_numbp50_final_div.tab", 'r')
counter = 0
seed = 0
gene_num = 1
d_col = {}
for aline in f_exons:
	aline1 = aline.strip('\n')
	aline2 = aline1.split('\t')
	#print (aline2)
	if aline2[0] == "exon":
		col = 0
		for x in aline2:
			d_col[x] = col
			col += 1
	else:
		gene_name = aline2[d_col["gene"]]
		exon_len = aline2[d_col["exon_len"]]
		inter_len = int(aline2[d_col["numbp50"]])*4
		rec_rate = str(round(float(aline2[d_col["avg_rec_rate"]]),3))
		repID = 1
		while repID <= num_rep:
			seed = seed + 1
			#check if this sim has already been done:
			try:
				f_check = open("/home/pjohri/DEMO_DFE_SINGEXON_HUMAN_V2/sim" + str(simID) + "/sim" + str(simID) + "_gene" + str(gene_num) + "_rep" + str(repID) + ".ms", 'r')
				f_check.close()
			except:
				counter = counter + 1
				try:
					f_csv.write(str(simID) + "," + str(gene_num) + "," + str(repID) + "," + str(seed) + "," + str(scaling_factor) + "," + str(Na) + "," + str(g_factor) + "," + str(f0) + "," + str(f1) + "," + str(f2) + "," + str(f3) + "," + str(rec_rate) + "," + str(inter_len) + "," + str(exon_len) + '\n')
				except:
					f_csv = open("/home/pjohri/SUBMITFILES/sim" + str(simID) + ".csv", 'w+')
					f_csv.write(str(simID) + "," + str(gene_num) + "," + str(repID) + "," + str(seed) + "," + str(scaling_factor) + "," + str(Na) + "," + str(g_factor) + "," + str(f0) + "," + str(f1) + "," + str(f2) + "," + str(f3) + "," + str(rec_rate) + "," + str(inter_len) + "," + str(exon_len) + '\n')
				f_script = open("/home/pjohri/demo_dfe_SingExon_human.slim", 'r')
				result = open("/home/pjohri/SLIMSCRIPTS/script_sim" + str(simID) + ".slim", 'w+')
				for line in f_script:
					line1 = line.replace("gen_burnin", str((10*Na) + (4*Na)))
					line2 = line1.replace("gen_stop", str((10*Na) + (4*Na) + 200))
					result.write(line2)
				#print("end of reading")
				result.close()
				f_script.close()
			repID = repID + 1
		gene_num = gene_num + 1

f_csv.close()

#copy and write the submit file with that simulation ID if there are any simulations left to perform in that simID:
if counter > 0:
	f_submit = open("/home/pjohri/SingExon_slim_queue_human_v2.submit", 'r')
	result_submit = open("/home/pjohri/SUBMITFILES/SingExon_slim_sim" + str(simID) + ".submit", 'w+')
	for line in f_submit:
        	line1 = line.replace("simID", "sim" + str(simID))
        	result_submit.write(line1)
	result_submit.close()
	f_submit.close()


print ("done")


