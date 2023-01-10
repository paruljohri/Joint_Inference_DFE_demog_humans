#This is to add the number of bases required for 50% and 75% recovery for a given recombination rate for specific exopns:

import sys
import math
import numpy as np

#Define variables and constants:                                                
g = 0.0 #rate of gene conversion                                                
u = 1.25*1e-8 #(*Mutation rate*)                                                                     
Ne = 20000.0 #(*Effective population size, only required to calculate expected nucleotide diversity under neutrality*)
pi = 4*Ne*u #(*Expected nucleotide diversity under neutrality*)                 
f0 = 0.22 #(*Proportion of effectively neutral mutations with 0 <= |2Nes| < 1 *) 
f1 = 0.27 #(*Proportion of weakly deleterious mutations with 1 <= |2Nes| < 10 *) 
f2 = 0.13 #(*Proportion of moderately deleterious mutations with 10 <= |2Nes| < 100 *)
f3 = 0.38 #(*Proportion of strongly deleterious mutations with |2Nes| >= 100 *)  
#(*Note that the number of classes can easily be increased to whatever is required to approximate the continuous DFE *)
h = 0.5 #(* dominance coefficient *)
t0 = 0.0                                                                        
t1 = h*(1/(2*Ne))                                                               
t1half = h*(5/(2*Ne)) #(* This is the cut-off value of 2Nes=5. This derivation assumes that all mutations with 2Nes<5 will not contribute to BGS *)
t2 = h*(10/(2*Ne))                                                              
t3 = h*(100/(2*Ne))                                                             
t4 = h*1.0

#defining functions required:
def calculate_exponent(t_start, t_end, posn, l, r):                                   
    U = u*l
    a = g + r*posn                                                              
    b = g + r*(posn + l)                                                        
    E1 = (U/(r*l*(1-a))) * (1 + (a/((1-a)*(t_end-t_start))) * math.log((a+(t_start*(1-a)))/(a + (t_end*(1-a)))))
    E2 = -1.0*(U/(r*l*(1-b)))*(1 + (b/((1-b)*(t_end-t_start)))*math.log((b + ((1-b)*t_start))/(b + ((1-b)*t_end))))
    E = E1 + E2                                                                 
    return (E)                                                                  
                                                                                
def calculate_pi(posn, l, r):                                                         
    E_f1 = calculate_exponent(t1half, t2, posn, l, r)                                 
    E_f2 = calculate_exponent(t2, t3, posn, l, r)                                     
    E_f3 = calculate_exponent(t3, t4, posn, l, r)                                     
    pi_posn = f0*pi + f1*0.5*pi + f1*0.5*pi*math.exp(-1.0*E_f1) + f2*pi*math.exp(-1.0*E_f2) + f3*pi*math.exp(-1.0*E_f3)
    return (pi_posn)

def fit_log_curve(l_posns, l_pi):                                               
    x_data = np.array(l_posns)                                                  
    y_data = np.array(l_pi)                                                     
    log_x_data = np.log(x_data)                                                 
    curve_fit = np.polyfit(log_x_data, y_data, 1)                               
    return (curve_fit[0], curve_fit[1])

def calculate_pi_segment(posn_start, posn_end, l, r):
    l_posns, l_pi = [], []                                                          
    j=int(posn_start)                                                               
    while j <= int(posn_end):                                                           
        l_posns.append(j)                                                           
        l_pi.append(calculate_pi(j, l, r))                                                
        j += 1                                                                      
    return(l_posns, l_pi)

d_col = {}
f_exon = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon1-6kb_CNE500_rec.tab", 'r')
result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon1-6kb_CNE500_rec_numbp50.tab", 'w+')
for line in f_exon:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] == "exon":
        result.write(line1 + '\t' + "numbp50" + '\t' + "numbp75" + '\n')
        col = 0
        for x in line2:
            d_col[x] = col
            col = col + 1
    else:
        print (line2[0])
        s_len = float(line2[d_col["exon_len"]])
        rec_check = line2[d_col["avg_rec_rate"]]
        if rec_check == "NA":
            result.write (line1 + '\t' + "NA" + '\t' + "NA" + '\n')
        else:
            rec_rate = float(line2[d_col["avg_rec_rate"]])*1e-8
            pi_segment = calculate_pi_segment(1, 100000, s_len, rec_rate)
            log_fit_parameters = fit_log_curve(pi_segment[0], pi_segment[1])
            slope = log_fit_parameters[0]
            intercept = log_fit_parameters[1]
            print (str(slope) + '\t' + str(intercept))
            pi_50 = float(calculate_pi(1, s_len, rec_rate)) + ((pi - float(calculate_pi(1, s_len, rec_rate)))*0.5)
            pi_75 = float(calculate_pi(1, s_len, rec_rate)) + ((pi - float(calculate_pi(1, s_len, rec_rate)))*0.75)
            result.write(line1 + '\t')
            try:    
                numbp50 = math.exp((pi_50-intercept)/slope)
                result.write(str(int(numbp50)) + '\t')
            except OverflowError:
                result.write("NA" + '\t')
            try:
                numbp75 = math.exp((pi_75-intercept)/slope)
                result.write(str(int(numbp75)) + '\n')
            except OverflowError:
                result.write("NA" + '\n')
f_exon.close()
result.close()
print ("done")




