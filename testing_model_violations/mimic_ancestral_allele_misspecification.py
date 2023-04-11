#This is to take 1% of all sites and change the singletons into 99% allel frequency:
# -inFolder -outFolder -simID -numGenes -numRep
import sys
import argparse
import random

#parsing user given constants                                                   
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-inFolder', dest = 'inFolder', action='store', nargs = 1, type = str, help = 'the name of the folder or simulation to run')
parser.add_argument('-outFolder', dest = 'outFolder', action='store', nargs = 1, type = str, help = 'the name of the folder in which output is placed')
#parser.add_argument('-interRegion', dest = 'interRegion', action='store', nargs = 1, type = str, help = '5p/3p; the intergenic region on which filtering is based')
parser.add_argument('-simID', dest = 'simulationID', action='store', nargs = 1, type = str, help = 'the name of the subfolder or simulation to run')
parser.add_argument('-numGenes', dest = 'numGenes', action='store', nargs = 1, type = int, help = 'number of genes')
parser.add_argument('-numRep', dest = 'numRep', action='store', nargs = 1, type = int, choices = range(1,10001), help = 'number of replicates of each gene')
args = parser.parse_args()                                                      
in_folder = args.inFolder[0]
out_folder = args.outFolder[0] 
simID = args.simulationID[0]                                              
num_genes = int(args.numGenes[0])
num_reps = args.numRep[0]

def convert_derived_to_ancestral(s_alleles):
    s_alleles_converted = ""
    for x in s_alleles:
        if x == '0':
            y = '1'
        elif x == '1':
            y = '0'
        else:
            print ("something weird about your genotypes")
        s_alleles_converted = s_alleles_converted + str(y)
    return (s_alleles_converted)

#read ms file with positions in their decimal forms                   
def read_subset_ms(f_ms, start, end):                               
    d_pos = {} #num of SNP -> scaled integer postion                            
    d_geno = {}#num of SNP -> genotypes                                         
    l_num_snps = []#list of number of SNPs                                      
    for line in f_ms:                                                           
        line1 = line.strip('\n')                                                
        if "positions" in line1:                                                
            line2 = line1.split()                                               
            i = 0                                                               
            for x in line2:                                                     
                if "position" not in x:                                         
                    if (float(x) >= float(start)) and (float(x) <= float(end)): 
                        #x_int = round(float(x)*(len_region-1))                  
                        d_pos[i] = x        
                        d_geno[i] = ""                                          
                        l_num_snps.append(i)                                    
                    i = i + 1                                                   
            #print(l_int_posn)                                                  
        elif "//" not in line and "segsites" not in line:                       
            for j in l_num_snps:
                d_geno[j] = d_geno[j] + line1[j]                                
    return (l_num_snps, d_pos, d_geno)

#write ms file when positions are only known as decimals:                      
def write_ms_file(RESULT, l_NUM_SNPS, d_POS, d_GENO):                
    result.write("//" + '\n')                                                   
    result.write("segsites: " + str(len(l_NUM_SNPS)) + '\n')                    
    result.write("positions:")
    for snp in l_NUM_SNPS:                                                      
        #POSN_NORM = round(float(d_POSN[snp])/float(len_region-1),7)             
        result.write(" " + str(d_POS[snp]))                                      
    result.write('\n')                                                          
    if len(l_NUM_SNPS) > 0:                                                     
        num_indv = len(d_GENO[l_NUM_SNPS[0]])                                     
        indv = 0                                                                
        while indv < num_indv:                                                  
            for snp in l_NUM_SNPS:                                              
                #print(indv)                                                    
                result.write(d_GENO[snp][indv])                                   
            result.write('\n')                                                  
            indv += 1                                                           
    return ("written to file")

#go through all simulation replicates
geneID = 0
num_singletons_modified = 0
while geneID < num_genes:
    geneID += 1
    repID=1
    while repID <= num_reps:
        #read ms file
        f_ms = open("/scratch/pjohri1/BgsDfeDemo_Human/" + in_folder + "/sim" + str(simID) + "/sim" + str(simID) + "_gene" + str(geneID) + "_rep" + str(repID) + ".ms", 'r')
        t_ms = read_subset_ms(f_ms, 0.0, 1.0)  
        f_ms.close() 
        d_gt = t_ms[2] #SNP num -> genotypes
        d_gt_new = {}
        for snp in d_gt.keys():
            #check if this snp is a dervied singleton:
            if d_gt[snp].count("1") == 1:
                coin_flip = random.uniform(0.0,1.0)
                if coin_flip <= 0.01:
                    num_singletons_modified += 1
                    new_gt = convert_derived_to_ancestral(d_gt[snp])
                    d_gt_new[snp] = new_gt
                else:
                    d_gt_new[snp] = d_gt[snp]
            else:
                d_gt_new[snp] = d_gt[snp]
        #write ms file
        print ("number of converted singletons: " + str(num_singletons_modified))
        result = open("/scratch/pjohri1/BgsDfeDemo_Human/" + out_folder + "/sim" + str(simID) + "_aa/sim" + str(simID) + "_gene" + str(geneID) + "_rep" + str(repID) + ".ms", 'w+')
        write_ms_file(result, t_ms[0], t_ms[1], d_gt_new)
        result.close()
        repID += 1

print ("done")
