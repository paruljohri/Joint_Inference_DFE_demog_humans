#This is to perform 1 replicate of each exon:
#Uses mutation rate for each exon separately
#This program will create a temporary slim script and save it's output, and then delete the script.
#python write_slim_scripts_SingExon_human.py -numRep 10 -simID 1
#go over starting ID to end ID. Check if a set has already been simulated. Then submit it if it hasn't been done.

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

#make a directory for results of this simualtion:
#os.system("mkdir /scratch/pjohri11/demo_dfe_SingExon_human/sim" + str(simID))

def average_of_list(l_values):
	s_sum, s_num = 0.0, 0.0
	for x in l_values:
		if x != "NA":
			s_sum = s_sum + float(x)
			s_num = s_num + 1.0
	return(s_sum/s_num)

#get parameters and ID of that simulation:
f_para = open("/home/pjohri1/BgsDfeDemo_Human/programs_model_violations/model_violations_parameters.txt", 'r')
d_col = {}
for line in f_para:
	line1 = line.strip('\n')
	line2 = line1.split()
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
		script_name = line2[d_col["script"]]
f_para.close()

#get mutation rate for every gene:
f_mutn = open("/home/pjohri1/BgsDfeDemo_Human/Humans/MutationRateMap/mutn_rate_single_exons.txt",'r')
d_mutn = {}
for line in f_mutn:
	line1 = line.strip('\n')
	line2 = line1.split('\t')
	if line2[0] != "exon":
		if line2[2] == "NA":
			s_rate = 1.25e-8
		else:
			s_rate = float(line2[2])*1.25
		try:
			d_mutn[line2[0]].append(s_rate)
		except:
			d_mutn[line2[0]] = []
			d_mutn[line2[0]].append(s_rate)
f_mutn.close()

#re-write the slim script only to change the burn-in period (that is generations)
f_script = open("/home/pjohri1/BgsDfeDemo_Human/programs_model_violations/" + script_name, 'r')
result = open("/scratch/pjohri1/BgsDfeDemo_Human/SLIMFILES/script_sim" + str(simID) + ".slim", 'w+')
for line in f_script:
	line1 = line.replace("gen_burnin", str((10*Na) + (4*Na)))
	line2 = line1.replace("gen_stop", str((10*Na) + (4*Na) + 200))
	result.write(line2)
result.close()
f_script.close()

#Write a bash script that functions as a submit script- it submits the geneIDs required to run:
f_submit = open("/scratch/pjohri1/BgsDfeDemo_Human/BASHFILES/sim" + str(simID) + "_submit.sh", 'w+')
f_submit.write("#!/bin/bash" + '\n')

#write a bash file for every simID+geneID with the full slim command line
f_exons = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon2-6kb_CNE500_rec_numbp50_final_div.tab", 'r')
counter = 0
seed = 0
gene_num = 1
d_col = {}
exon_num = 0
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
		exon_num += 1
		gene_name = aline2[d_col["gene"]]
		s_exon = "exon" + str(exon_num) + "_" + gene_name
		exon_len = aline2[d_col["exon_len"]]
		inter_len = int(aline2[d_col["numbp50"]])*4
		rec_rate = str(round(float(aline2[d_col["avg_rec_rate"]]),3))
		#check if any of the 10 replicates of this gene are still left to be simulated:
		check = 0
		repID = 1
		while repID <= num_rep:
			try:
				f_check = open("/scratch/pjohri1/BgsDfeDemo_Human/ModelViolations/sim" + str(simID) + "/sim" + str(simID) + "_gene" + str(gene_num) + "_rep" + str(repID) + ".ms", 'r')
				f_check.close()
			except:
				check = 1
				counter += 1 #counts total number of reps in this simID that are still left to be simulated
			repID = repID + 1
		if check == 1:#at least one replicate is missing
			f_bash = open("/scratch/pjohri1/BgsDfeDemo_Human/BASHFILES/sim" + str(simID) + "_gene" + str(gene_num) + ".sh", 'w+')
			f_bash.write("#!/bin/bash" + '\n')
			f_bash.write("repID=$1" + '\n')
			f_bash.write("seed=$2" + '\n')
			f_bash.write("slim -d d_seed=${seed} -d d_Nanc_scaled=" + str(Na) + " -d d_scaling_factor=" + str(scaling_factor) + " -d d_growth=" + str(g_factor) + " -d d_f0=" + str(f0) + " -d d_f1=" + str(f1) + " -d d_f2=" + str(f2) + " -d d_f3=" + str(f3) + " -d d_rec_rate=" + str(rec_rate) + " -d d_mut_rate=" + str(average_of_list(d_mutn[s_exon])) + " -d d_inter_len=" + str(inter_len) + " -d d_exon_len=" + str(exon_len) + " -d \"d_simID='" + str(simID) + "'\" -d \"d_geneID='" + str(gene_num) + "'\" -d \"d_repID='${repID}'\" -d \"d_folder='/scratch/pjohri1/BgsDfeDemo_Human/ModelViolations/sim" + str(simID) + "'\" /scratch/pjohri1/BgsDfeDemo_Human/SLIMFILES/script_sim" + str(simID) + ".slim" + '\n')
			f_submit.write("sbatch /home/pjohri1/BgsDfeDemo_Human/programs_model_violations/run_SingExon_slim_array_human_agave.sh " + str(simID) + " " + str(gene_num) + '\n')#writing this geneID in the submit file
			repID = repID + 1
		gene_num = gene_num + 1
f_submit.close()

#if all reps are done then delete the slim script created for that simID:
if counter == 0:
	os.system("rm /scratch/pjohri1/BgsDfeDemo_Human/SLIMFILES/script_sim" + str(simID) + ".slim")


print ("done")


