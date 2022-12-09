#This is to filter VCFs by the following criteria:
#Exclude any variants not present in the accessibility regions.
#Exlcude variants present in phastcons elements.
#Note down these sites in a file so they can be used later to filter variants and divergence and also use as a denominator:
#Changing this to be per chromosome.

import sys

s_chrom = sys.argv[1]
l_windows = ["linked1", "linked2", "linked3", "linked4"]

#get chromosome names                                                           
d_chrom_num = {}                                                                    
f_info = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_info.txt", 'r')
for line in f_info:                                                             
    line1 = line.strip('\n')                                                    
    line2 = line1.split('\t')                                                   
    d_chrom_num[line2[6]] = line2[9]                                                
f_info.close()

#get the coordinates of the exons:
print("reading exon coordinates")
f_exons = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon2-6kb_CNE500_rec_numbp50_final_div.tab", 'r')
d_col = {}                                                                      
l_exons = []                                                                    
exon_num = 1
d_chr = {}                                                                    
d_strand, d_numbp50 = {}, {}                                                    
#d_exon_start, d_exon_end = {}, {}
d_coord_start = {}
d_coord_end = {}
l_regions = ["exon", "5p", "3p"]
for region in l_regions:
    d_coord_start[region] = {}
    d_coord_end[region] = {}
for line in f_exons:                                                            
    lineA = line.strip('\n')                                                    
    lineB = lineA.split('\t')                                                   
    if lineB[0] == "exon":                                                      
        col_num = 0                                                             
        for x in lineB:                                                         
            d_col[x] = col_num                                                  
            col_num += 1
    else:
        s_gene = lineB[d_col["gene"]]                                           
        s_exon = "exon" + str(exon_num) + "_" + s_gene
        s_chr = d_chrom_num[lineB[d_col["chrom"]]]
        if s_chrom == s_chr:
            l_exons.append(s_exon)
            d_chr[s_exon] = s_chr
            d_strand[s_exon] = lineB[d_col["strand"]]                               
            d_coord_start["exon"][s_exon] = int(lineB[d_col["exon_start"]])
            d_coord_end["exon"][s_exon] = int(lineB[d_col["exon_end"]])
            d_numbp50[s_exon] = int(lineB[d_col["numbp50"]])
        exon_num += 1                                                           
f_exons.close()

#This is just for troubleshooting. Comment it out later!
#l_exons = ["exon9_ANKRD13C"]


def get_window(SITE, START, END, NUMBP50):
    if int(START) <= int(SITE) <= int(START)+int(NUMBP50)-1:
        WIN = "linked1"
    elif int(START)+int(NUMBP50) <= int(SITE) <= int(START)+(2*int(NUMBP50))-1:
        WIN = "linked2"
    elif int(START)+(2*int(NUMBP50)) <= int(SITE) <= int(START)+(3*int(NUMBP50))-1:
        WIN = "linked3"
    elif int(START)+(3*int(NUMBP50)) <= int(SITE) <= int(END):
        WIN = "linked4"
    else:
        print ("error in finding the correct window.")
    return(WIN)

d_windows = {}
def assign_windows(START, END, NUMBP50, DIR):   
    SITE = int(START)
    while SITE <= int(END):                                
        if int(START) <= int(SITE) <= int(START)+int(NUMBP50)-1:
            if DIR == "exon_at_start":
                WIN = "linked1"
            elif DIR == "exon_at_end":
                WIN = "linked4"                                                         
        elif int(START)+int(NUMBP50) <= int(SITE) <= int(START)+(2*int(NUMBP50))-1: 
            if DIR == "exon_at_start":
                WIN = "linked2"
            elif DIR == "exon_at_end":
                WIN = "linked3"
        elif int(START)+(2*int(NUMBP50)) <= int(SITE) <= int(START)+(3*int(NUMBP50))-1:
            if DIR == "exon_at_start":
                WIN = "linked3"
            elif DIR == "exon_at_end":
                WIN = "linked2"
        elif int(START)+(3*int(NUMBP50)) <= int(SITE) <= int(END):
            if DIR == "exon_at_start":
                WIN = "linked4"
            elif DIR == "exon_at_end":
                WIN = "linked1"
        else:                                                                       
            print ("error in finding the correct window.")
        d_windows[SITE] = WIN
        SITE += 1
    return()


#Get coordinates of different regions:
for exon in l_exons:
    if d_strand[exon] == "+":
        d_coord_start["5p"][exon] = int(d_coord_start["exon"][exon]) - 4.0*int(d_numbp50[exon])
        d_coord_end["5p"][exon] = int(d_coord_start["exon"][exon]) - 1
        assign_windows(d_coord_start["5p"][exon], d_coord_end["5p"][exon], d_numbp50[exon], "exon_at_end")
        d_coord_start["3p"][exon] = int(d_coord_end["exon"][exon]) + 1
        d_coord_end["3p"][exon] = int(d_coord_end["exon"][exon]) + 4.0*int(d_numbp50[exon])
        assign_windows(d_coord_start["3p"][exon], d_coord_end["3p"][exon], d_numbp50[exon], "exon_at_start")
    elif d_strand[exon] == "-":
        d_coord_start["5p"][exon] = int(d_coord_end["exon"][exon]) + 1
        d_coord_end["5p"][exon] = int(d_coord_end["exon"][exon]) + 4.0*int(d_numbp50[exon])
        assign_windows(d_coord_start["5p"][exon], d_coord_end["5p"][exon], d_numbp50[exon], "exon_at_start")
        d_coord_start["3p"][exon] = int(d_coord_start["exon"][exon]) - 4.0*int(d_numbp50[exon])
        d_coord_end["3p"][exon] = int(d_coord_start["exon"][exon]) - 1
        assign_windows(d_coord_start["3p"][exon], d_coord_end["3p"][exon], d_numbp50[exon], "exon_at_end")
#print(d_coord_start)
#print(d_coord_end)

#get the starting set of sites for each exon and each region:
d_sites = {}
d_present_sites = {}
d_acc_sites, d_cons_sites = {}, {}
for region in l_regions:                                                        
    d_sites[region] = {}
    d_acc_sites[region] = {}
    d_cons_sites[region] = {}
    for exon in l_exons:
        i = d_coord_start[region][exon]
        while i <= d_coord_end[region][exon]:
            d_present_sites[i] = "True"
            try:
                d_sites[region][exon].append(i)
            except:
                d_sites[region][exon] = []
                d_sites[region][exon].append(i)
            i = i + 1

#read the coordinates of the accessible parts:
print("adding accessible sites")
f_acc = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Accessibility/strictAccs_1000G_all_Sep2022", 'r')
for line in f_acc:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line[0] != "#":
        s_chr = line2[0]
        if s_chr == s_chrom:
            s_begin = int(line2[1])
            s_end = int(line2[2])
            for site in range(s_begin, s_end+1):
                if d_present_sites.get(int(site), "NA") == "True":
                    d_acc_sites[int(site)] = "A"
#print("all sites in accessible regions")
#print (d_sites)

#subtract the sites in the phastcons elements
print("subtracting phastcons elements")
f_cons = open("/home/pjohri1/BgsDfeDemo_Human/Humans/phastCons/phastcons_humans_100vertebrates.txt", 'r')
for line in f_cons:
    line1 = line.strip('\n')
    line2 = line.split('\t')
    if line2[0] != "bin":
        s_chr = line2[1]
        if s_chr == s_chrom:
            s_begin = int(line2[2])+1 #position is 0-based
            s_end = int(line2[3])+1 #position is 0-based
            for site in range(s_begin, s_end+1):
                if d_present_sites.get(int(site), "NA") == "True":
                    d_cons_sites[int(site)] = "C"

#Write out the sites:
if s_chrom == "chr1":
    result_num_sites = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/FILTERED_SITES/num_of_filtered_sites.txt", 'w+')
    result_num_sites.write("exon" + '\t' + "region" + '\t' + "window" + '\t' + "num_tot_sites" + '\t' + "num_filtered_sites" + '\n')
else:
    result_num_sites = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/FILTERED_SITES/num_of_filtered_sites.txt", 'a+')
print("writing files")
for exon in l_exons:
    for region in l_regions:
        s_exon_tot = 0
        result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/FILTERED_SITES/" + exon + "_" + region + ".sites", 'w+')
        d_num_sites = {}
        for win in l_windows:
            d_num_sites[win] = 0
        for site in d_sites[region][exon]:
            if d_acc_sites.get(int(site), "NA") == "A": #site is accessible
                if region == "exon":#if in coding region, don't eliinate phastcons elements
                    result.write(str(round(site)) + '\n')
                    s_exon_tot += 1
                else:#if region in not coding, then check for phastcons elements:
                    if d_cons_sites.get(int(site), "NA") == "NA": #does not belong to a phastcons element
                        result.write(str(round(site)) + '\n')
                    d_num_sites[d_windows[site]] += 1
        result.close()
        #write out the number of sites:
        if region == "exon":
            result_num_sites.write(exon + '\t' + region + '\t' + "func"  + '\t' + str(s_exon_tot) + '\n')
        else:
            for win in l_windows:
                result_num_sites.write(exon + '\t' + region + '\t' + win  + '\t' + str(d_num_sites[win]) + '\n')
result_num_sites.close()
print ("done")



