#This is to find folders that have all simulations complete and then zip them:
#This script is best run directly in the folder.

import sys
import os

s_start = sys.argv[1]
s_end = sys.argv[2]

def get_num_lines_in_file(f_name):
	num_lines = 0
	for line in f_name:
		num_lines += 1
	return(num_lines)

simID = int(s_start)
while simID <= int(s_end):
	print (simID)
	subfolder = "sim" + str(simID)
	os.system("ls " + subfolder + "/*.ms > tmp.txt")
	f_list = open("tmp.txt", 'r')
	numLines = get_num_lines_in_file(f_list)
	f_list.close()
	if numLines == 465:
		print("zipping folder")
		os.system("zip -r " + subfolder + ".zip " + subfolder)
		if os.path.getsize(subfolder + ".zip") > 0:#OR os.stat("filename").st_size > 0
			os.system("rm -r " + subfolder)
			#also delete the log files:
			print("deleting the log files")
			os.system("rm /home/pjohri/LOGS/slimtest_sim" + str(simID) + "_gene*")
	simID += 1
print ("done")

