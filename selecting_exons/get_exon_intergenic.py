#This is to get coordinates of first and last exon, the upstream and downstream intergenic region, also get corresponding full gene length and % of introns in it:
#Getting the longest UTR regions as coordinates
#Currently this program ignores all genes at the end of the chromosomes. Include that if needed.
#When there are overlapping exonso or genes, it takes the longest exon size.

import sys

def get_geneID(l_line2):
    for x in l_line2[8].split(";"):
        if "gene=" in x:
            geneID0 = x
    geneID = geneID0.replace("gene=", "")
    return geneID
def get_exon_type(l_line2):
    s_type = "NA"
    for x in l_line2[8].split(";"):
        if "gbkey" in x:
            s_type = x.replace("gbkey=", "")
    return s_type
def get_gene_type(l_line2):                                                     
    s_type = "NA"                                                               
    for x in l_line2[8].split(";"):                                             
        if "gene_biotype" in x:                                                        
            s_type = x.replace("gene_biotype=", "")                                    
    return s_type

#get gene type for each gene:
f_ann = open("/home/pjohri1/BgsDfeDemo_Human/Humans/RefGenome_hg19/GRCh37_latest_genomic_sorted.gff", 'r')
d_exon_type = {}
d_gene_type = {}
for line in f_ann:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[2] == "exon":
        tmp1 = get_geneID(line2).split(":")
        tmp1.pop()
        geneID = ":".join(tmp1)
        d_exon_type[geneID] = get_exon_type(line2)
        if geneID == "snoRNA:kis-a":
            print (line2)
    elif line2[2] == "gene":
        geneID = get_geneID(line2)
        d_gene_type[geneID] = get_gene_type(line2)
f_ann.close()

#Now re-read it:        
f_ann = open("/home/pjohri1/BgsDfeDemo_Human/Humans/RefGenome_hg19/GRCh37_latest_genomic_sorted.gff", 'r')
result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_intergenic.tab", 'w+')
l_exons = []
d_exon_start, d_exon_end = {}, {}
d_chrom = {}
d_chrom[""] = ""
d_strand = {}
d_strand[""] = ""
d_inter_5p_start, d_inter_3p_start = {}, {}
d_inter_5p_end, d_inter_3p_end = {}, {}
exonID0 = ""
geneID0 = ""
exon_start = ""
exon_end = "1"
d_exon_num = {}
d_len_exon = {}
d_gene = {}
for line in f_ann:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[2]=="exon":#use this to define intergenic
        geneID = get_geneID(line2)
        if d_exon_num.get(geneID, "NA") == "NA":
            d_exon_num[geneID] = 1
        else:
            d_exon_num[geneID] += 1
        #exonID = geneID + "_" + str(d_exon_num[geneID])
        exonID = line2[0] + "-" + line2[3] + "-" + line2[4]
        if exonID != exonID0:
            s_len = int(line2[4]) - int(line2[3]) + 1
            if s_len >= 50:#any element larger than 50 bp shoudl be considered as an exon
                #print (exonID)
                d_gene[exonID] = geneID
                d_len_exon[exonID] = int(line2[4]) - int(line2[3])
                d_chrom[exonID] = line2[0]
                d_strand[exonID] = line2[6]
                d_exon_start[exonID] = line2[3]
                d_exon_end[exonID] = line2[4]
                l_exons.append(exonID)
                if d_strand[exonID] == "+":
                    if d_chrom[exonID] == d_chrom[exonID0]:
                        d_inter_5p_start[exonID] = int(exon_end) + 1
                    else:
                        d_inter_5p_start[exonID] = "begin"
                    d_inter_5p_end[exonID] = int(line2[3]) - 1
                elif d_strand[exonID] == "-":
                    if d_chrom[exonID] == d_chrom[exonID0]:
                        d_inter_3p_start[exonID] = int(exon_end) + 1
                    else:
                        d_inter_3p_start[exonID] = "begin"
                    d_inter_3p_end[exonID] = int(line2[3]) - 1
                if d_strand[exonID0] == "+":
                    d_inter_3p_start[exonID0] = int(exon_end) + 1
                    if d_chrom[exonID] == d_chrom[exonID0]:
                        d_inter_3p_end[exonID0] = int(line2[3]) - 1
                    else:
                        d_inter_3p_end[exonID0] = "end"
                elif d_strand[exonID0] == "-":
                    d_inter_5p_start[exonID0] = int(exon_end) + 1
                    if d_chrom[exonID] == d_chrom[exonID0]:
                        d_inter_5p_end[exonID0] = int(line2[3]) - 1
                    else:
                        d_inter_5p_end[exonID0] = "end"
                if exonID0 == "NC_000001.10-170120519-170120603":
                    print (d_chrom[exonID0])
                    print (d_strand[exonID0])
                    print (d_inter_5p_start[exonID0])
                    print (d_inter_5p_end[exonID0])
                    print (d_inter_3p_start[exonID0])
                    print (d_inter_3p_end[exonID0])
                exon_start = line2[3]
                exon_end  = line2[4]
                exonID0 = exonID
#adding intergenic information to the last gene:
if d_strand[exonID] == "+":                                     
    if d_chrom[exonID] == d_chrom[exonID0]:                     
        d_inter_5p_start[exonID] = int(exon_end) + 1            
    else:                                                       
        d_inter_5p_start[exonID] = "begin"                      
    d_inter_5p_end[exonID] = int(line2[3]) - 1                  
elif d_strand[exonID] == "-":                                   
    if d_chrom[exonID] == d_chrom[exonID0]:                     
        d_inter_3p_start[exonID] = int(exon_end) + 1            
    else:                                                       
        d_inter_3p_start[exonID] = "begin"                      
    d_inter_3p_end[exonID] = int(line2[3]) - 1                  
if d_strand[exonID0] == "+":                                    
    d_inter_3p_start[exonID0] = int(exon_end) + 1               
    if d_chrom[exonID] == d_chrom[exonID0]:                     
        d_inter_3p_end[exonID0] = int(line2[3]) - 1             
    else:                                                       
        d_inter_3p_end[exonID0] = "end"                         
elif d_strand[exonID0] == "-":                                  
    d_inter_5p_start[exonID0] = int(exon_end) + 1               
    if d_chrom[exonID] == d_chrom[exonID0]:                     
        d_inter_5p_end[exonID0] = int(line2[3]) - 1             
    else:                                                       
        d_inter_5p_end[exonID0] = "end"

result.write("exon" + '\t' + "gene" + '\t' + "chrom" + '\t' + "exon_type" + '\t' + "strand" + '\t' + "exon_len" + '\t' + "exon_start"  + '\t' + "exon_end" + '\t' + "intergenic_5p_start" + '\t' + "intergenic_5p_end" + '\t' + "intergenic_3p_start" + '\t' + "intergenic_3p_end" + '\t' + "intergenic_5p_len" + '\t' + "intergenic_3p_len" + '\n')
l_default = ["NA"]
for exon in l_exons:
    print(exon)
    print (d_inter_5p_start[exon])
    print (d_inter_3p_start[exon])
    print (d_inter_5p_end[exon])
    print (d_inter_3p_end[exon])
    result.write(exon + '\t' + d_gene[exon] + '\t' + d_chrom[exon] + '\t' + d_gene_type[d_gene[exon]] + '\t' + d_strand[exon] + '\t' + str(d_len_exon[exon]) + '\t' + str(d_exon_start[exon])  + '\t' + str(d_exon_end[exon]) + '\t' + str(d_inter_5p_start[exon]) + '\t' + str(d_inter_5p_end[exon]) + '\t' + str(d_inter_3p_start[exon]) + '\t' + str(d_inter_3p_end[exon]))
    if isinstance(d_inter_5p_start[exon], str) == True or isinstance(d_inter_5p_end[exon], str) == True:
        result.write('\t' + "NA")
    else:
        result.write('\t' + str(d_inter_5p_end[exon]-d_inter_5p_start[exon]+1))
    if isinstance(d_inter_3p_start[exon], str) == True or isinstance(d_inter_3p_end[exon], str) == True:
        result.write('\t' + "NA" + '\n')
    else:
        result.write('\t' + str(d_inter_3p_end[exon]-d_inter_3p_start[exon]+1) + '\n')
    #print (d_exon_type[gene])
#   result.write(gene + '\t' + d_chrom[gene] + '\t' + d_exon_type.get(gene, "NA") + '\t' + d_strand[gene] + '\t' + str(d_len_gene[gene]))
    #Figuring out the first exon (this is being doen because the annotation file doesnt always have this in order)
#   s_start = ""
#   s_end = ""
    #check if exon coordinates are provided
#   if d_exon_posn.get(gene, "NA") == "NA" and d_cds_posn.get(gene, "NA") != "NA":
        #print (gene)
        #print (d_cds_posn[gene])
        #Finding first CDS:
#       for x in d_cds_posn[gene]:
#           s_start1 = x.split("-")[0]
#           if s_start == "":
#               s_start = s_start1[:]
#               s_end = x.split("-")[1]
#           else:
#               if int(s_start1) < int(s_start):
#                   s_start = s_start1[:]
#                   s_end = x.split("-")[1]
        #get utr information and add to CDS:
#       if d_strand[gene] == "+":
#           start_5p = str(min(d_u5_start.get(gene, [s_start])))
#           end_5p = s_end[:]#end of cds
#       elif d_strand[gene] == "-":
#           start_3p = str(min(d_u3_start.get(gene, [s_start])))
#           end_3p = s_end[:]
        #Finding last CDS:
#       s_start = ""
#       s_end = ""
#       if d_cds_posn.get(gene, "NA")=="NA":
#           print (gene + " check something weird")
#           s_start, s_end = "NA", "NA"
#       else:
#           for x in d_cds_posn[gene]:
#               s_end1 = x.split("-")[1]
#               if s_end == "":
#                   s_end = s_end1[:]
#                   s_start = x.split("-")[0]
#               else:
#                   if int(s_end1) > int(s_end):
#                       s_end = s_end1[:]
#                       s_start = x.split("-")[0]
#       if d_strand[gene] == "+":
#           start_3p = s_start[:]
#           end_3p = str(max(d_u3_end.get(gene, [s_end])))
#       elif d_strand[gene] == "-":
#           start_5p = s_start[:]
#           end_5p = str(max(d_u5_end.get(gene, [s_end])))
            
#   elif d_exon_posn.get(gene, "NA") != "NA":#using exon information
        #print (gene)
        #print (d_exon_posn[gene])
#       for x in d_exon_posn[gene]:
#           s_start1 = x.split("-")[0]
#           if s_start == "":
#               s_start = s_start1[:]
#               s_end = x.split("-")[1]
#           else:
#               if int(s_start1) < int(s_start):
#                   s_start = s_start1[:]
#                   s_end = x.split("-")[1]
#       if d_strand[gene] == "+":
            #print (s_start)
#           start_5p = s_start[:]
#           end_5p = s_end[:]
#       elif d_strand[gene] == "-":
#           start_3p = s_start[:]
#           end_3p = s_end[:]
        #Figuring out the last exon
#       s_start = ""
#       s_end = ""
#       if d_exon_posn.get(gene, "NA")=="NA":
#           print (gene + " check something weird")
#           s_start, s_end = "NA", "NA"
#       else:
#           for x in d_exon_posn[gene]:
#               s_end1 = x.split("-")[1]
#               if s_end == "":
#                   s_end = s_end1[:]
#                   s_start = x.split("-")[0]
#               else:
#                   if int(s_end1) > int(s_end):
#                       s_end = s_end1[:]
#                       s_start = x.split("-")[0]
#       if d_strand[gene] == "+":
#           start_3p = s_start[:]
#           end_3p = s_end[:]
#       elif d_strand[gene] == "-":
#           start_5p = s_start[:]
#           end_5p = s_end[:]
#   result.write('\t' + start_5p + '\t' + end_5p + '\t' + start_3p + '\t' + end_3p)
#   result.write('\t' + ",".join(d_intron_posn.get(gene, l_default)))
#   l_intron_len = []
#   for x in d_intron_posn.get(gene, l_default):
#       if x == "NA":
#           l_intron_len = ["NA"] 
#       else:
#           s_len = int(x.split("-")[1]) - int(x.split("-")[0]) + 1
#           l_intron_len.append(str(s_len))
#   result.write('\t' + ",".join(l_intron_len) + '\t')
#   if d_strand[gene] == "+":
#       result.write(str(l_intron_len[0]) + '\t' + str(l_intron_len[len(l_intron_len)-1]) + '\t')
#   elif d_strand[gene] == "-":
#       result.write(str(l_intron_len[len(l_intron_len)-1]) + '\t' + str(l_intron_len[0]) + '\t')
#   if d_exon_type.get(gene, "NA") == "mRNA":
#       result.write(d_inter_5p.get(gene, "NA") + '\t' + d_inter_3p.get(gene, "NA") + '\t')
#       if "begin" in d_inter_5p[gene] or "begin" in d_inter_3p[gene] or "end" in d_inter_5p[gene] or "end" in d_inter_3p[gene]:
#           result.write("NA" + '\t' + "NA" + '\n')
#       else:
#           result.write(str(int(d_inter_5p[gene].split("-")[1]) - int(d_inter_5p[gene].split("-")[0]) + 1) + '\t' + str(int(d_inter_3p[gene].split("-")[1]) - int(d_inter_3p[gene].split("-")[0]) + 1) + '\n') 
#   else:
#       result.write("NA" + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\n')
result.close()
f_ann.close()

print ("done")

