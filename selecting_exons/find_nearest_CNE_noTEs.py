#This is to add other elements in the intergenic between exons:
#Only use elements that are larger than 50 bp unless they are TEs.

import sys
chromname_ucsc = sys.argv[1]
element_len_cutoff = 500

#Reading conversion of chromosome names:
print ("reading chromosome info")
d_ncbi_chromname, d_ucsc_chromname = {}, {}
d_elements = {}
f_info = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_info.txt", 'r')
d_col = {}
for line in f_info:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] == "Sequence-Name":
        col = 0
        for x in line2:
            d_col[x] = col
            col = col + 1
    else:
        d_ncbi_chromname[line2[d_col["UCSC-style-name"]]] = line2[d_col["RefSeq-Accn"]]
        d_ucsc_chromname[line2[d_col["RefSeq-Accn"]]] = line2[d_col["UCSC-style-name"]]
        d_elements[line2[d_col["UCSC-style-name"]]] = []
f_info.close()
print (d_ncbi_chromname)

#Reading different annotation types:
d_chrom = {}
d_start = {}
d_end = {}

#Reading SNO and miRNAs:
print ("reading sno and miRNAs")
f_sno = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/sno_miRNA.txt", 'r')
d_col = {}
for line in f_sno:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] == "#bin":
        col = 0
        for x in line2:
            d_col[x] = col
            col = col + 1
    else:
        if line2[d_col["chrom"]] == chromname_ucsc:
            s_len = int(line2[d_col["chromEnd"]]) - int(line2[d_col["chromStart"]]) + 1
            if s_len >= element_len_cutoff:
                s = line2[d_col["name"]] + "-" + line2[d_col["chromStart"]] + "-" + line2[d_col["chromEnd"]]
                d_chrom[s] = line2[d_col["chrom"]]
                d_start[s] = int(line2[d_col["chromStart"]])
                d_end[s] = int(line2[d_col["chromEnd"]])
                d_elements[d_chrom[s]].append(s)
f_sno.close()

#phastcons elements:
print ("reading phastCons")
f_phast = open("/home/pjohri1/BgsDfeDemo_Human/Humans/phastCons/phastcons_humans_100vertebrates.txt", 'r')
d_col = {}
for line in f_phast:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] == "bin":
        col = 0
        for x in line2:
            d_col[x] = col
            col = col + 1
    else:
        if line2[d_col["chrom"]] == chromname_ucsc:
            s_len = int(line2[d_col["chromEnd"]]) - int(line2[d_col["chromStart"]]) + 1
            if s_len >= element_len_cutoff:
                s = "phastConsE" + "-" + line2[d_col["chromStart"]] + "-" + line2[d_col["chromEnd"]]
                d_chrom[s] = line2[d_col["chrom"]]
                d_start[s] = int(line2[d_col["chromStart"]])
                d_end[s] = int(line2[d_col["chromEnd"]])
                if d_chrom[s] in d_elements.keys():
                    d_elements[d_chrom[s]].append(s)
f_phast.close()
def get_start(val):
    return(d_start[val])
def get_end(val):
    return(d_end[val])
print ("going through intergenic table")
#Go through the table with intergenics between exons:
f_int = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_intergenic.tab", 'r')
result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_" + str(element_len_cutoff) + "bp/GRCh37_exons_firstelement_" + chromname_ucsc + "_" + str(element_len_cutoff) + "bp.tab", 'w+')
d_col = {}
for line in f_int:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] == "exon":
        result.write(line1 + '\t' + "5p_end" + '\t' + "5p_elements" + '\t' + "3p_end" + '\t' + "3p_elements" + '\n')
        col = 0
        for x in line2:
            d_col[x] = col
            col += 1
    else:
        print (line2[0])
        l_5p, l_3p = [], []
        s_chrom = line2[0].split("-")[0]
        if "begin" not in line1 and "end" not in line1:
            if d_ucsc_chromname.get(s_chrom, "NA") == chromname_ucsc and int(line2[d_col["intergenic_5p_len"]]) > 5000 and int(line2[d_col["intergenic_3p_len"]]) > 5000:
                result.write(line1)
                for element in d_elements[d_ucsc_chromname[s_chrom]]:
                    if line2[d_col["intergenic_5p_start"]] != "begin" and line2[d_col["intergenic_5p_end"]] != "end":
                        if d_start[element] >= int(line2[d_col["intergenic_5p_start"]]):
                            if d_end[element] <= int(line2[d_col["intergenic_5p_end"]]):
                                l_5p.append(element)
                    if line2[d_col["intergenic_3p_start"]] != "begin" and line2[d_col["intergenic_3p_end"]] != "end":
                        if d_start[element] >= int(line2[d_col["intergenic_3p_start"]]):
                            if d_end[element] <= int(line2[d_col["intergenic_3p_end"]]):
                                l_3p.append(element)
                if len(l_5p) == 0:
                    result.write('\t' + "NA" + '\t' + "none")
                else:
                    #Find the nearest element and write that down if the distnace between the nearest element depending on the strand
                    if line2[d_col["strand"]] == "+":
                        l_5p.sort(key=get_end)
                        elementID = l_5p.pop()
                        result.write('\t' + str(d_end[elementID]) + '\t' + elementID)
                    elif line2[d_col["strand"]] == "-":
                        l_5p.sort(key=get_start)
                        elementID = l_5p[0]
                        result.write('\t' + str(d_start[elementID]) + '\t' + elementID)
                if len(l_3p) == 0:
                    result.write('\t' + "NA" + '\t' + "none" + '\n')
                else:
                    if line2[d_col["strand"]] == "+":
                        l_3p.sort(key=get_start)
                        elementID = l_3p[0]
                        result.write('\t' + str(d_start[elementID]) + '\t' + elementID + '\n')
                    elif line2[d_col["strand"]] == "-":
                        l_3p.sort(key=get_end)
                        elementID = l_3p.pop()
                        result.write('\t' + str(d_end[elementID]) + '\t' + elementID + '\n')

result.close()
f_int.close()
print ("done")
















