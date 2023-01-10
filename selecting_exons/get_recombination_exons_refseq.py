#This is to get the recombination rate of exons:

import sys

def average(l_val):
	s_sum, s_tot = 0, 0
	for x in l_val:
		if x != "NA":
			s_tot = s_tot + 1
			s_sum = s_sum + float(x)
	if s_tot == 0:
		s_avg = "NA"
	else:
		s_avg = str(float(s_sum)/float(s_tot))
	return (s_avg)

#storing recombination rates:
d_rec = {}
d_chr = {}
d_start, d_end = {}, {}
f_rec = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Recombination/HapMap_24YRI_hg19", 'r')
for line in f_rec:
	line1 = line.strip('\n')
	if line1[0] != "#":
		line2 = line1.split('\t')
		s = line2[0] + "_" + line2[1] + "_" + line2[2]
		d_chr[s] = line2[0]
		d_start[s] = line2[1]
		d_end[s] = line2[2]
		d_rec[s] = line2[3]
f_rec.close()

#going throug exons:
f_exon = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon1-6kb_CNE500.tab", 'r')
result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon1-6kb_CNE500_rec.tab", 'w+')
d_col = {}
for line in f_exon:
	line1 = line.strip('\n')
	line2 = line1.split('\t')
	if line2[0] == "exon":
		result.write(line1 + '\t' + "rec_region_start" + '\t' + "rec_region_end" + '\t' + "avg_rec_rate" + '\n')
		col = 0
		for x in line2:
			d_col[x] = col
			col = col + 1
	else:
		s_start = int(line2[d_col["exon_start"]])
		s_end = int(line2[d_col["exon_end"]])
		s_chr = line2[d_col["chrom"]]
		gene = line2[d_col["gene"]]
		print (gene)
		#i = 0
		#while i < len(l_starts):
		#	 check = 0
			#print (l_starts[i] + '\t' + l_ends[i])
		#	 s_start = int(l_starts[i])
			#s_end = int(l_ends[i])
		l_rates = []
		l_coords = []
		for x in d_start.keys():
		#		 if d_chr[x] == s_chr:
			if int(d_start[x]) <= s_start:
				if s_end <= int(d_end[x]):#exon lies within a recomb region
					l_rates.append(d_rec[x])
					l_coords.append(int(d_start[x]))
					l_coords.append(int(d_end[x]))
			elif s_start <= int(d_start[x]):
				if int(d_end[x]) <= s_end:#recomb region lies within the exon
					l_rates.append(d_rec[x])
					l_coords.append(int(d_start[x]))
					l_coords.append(int(d_end[x]))
			elif int(d_start[x]) <= s_start:
				if s_start <= int(d_end[x]):#recomb region overlaps start of exon
					l_rates.append(d_rec[x])
					l_coords.append(int(d_start[x]))
					l_coords.append(int(d_end[x]))
			elif int(d_start[x]) <= s_end:
				if s_end <= int(d_end[x]):#recomb region overlaps end of exon
					l_rates.append(d_rec[x])
					l_coords.append(int(d_start[x]))
					l_coords.append(int(d_end[x]))
		l_coords.sort()
		if len(l_coords) >= 2:
			result.write(line1 + '\t' + str(l_coords[0]) + '\t' + str(l_coords.pop()) + '\t' + str(average(l_rates)) + '\n')
		else:
			result.write(line1 + '\t' + "NA" + '\t' + "NA" + '\t' + str(average(l_rates)) + '\n')

f_exon.close()
result.close()
print ("done")



