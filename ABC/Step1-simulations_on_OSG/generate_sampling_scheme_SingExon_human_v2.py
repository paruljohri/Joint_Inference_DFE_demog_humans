#This is to generate the sampling scheme for f0, f1, f2, f3, Nanc and Ncur:
#N is picked with loguniform distribution between upper and lower bound
#Scaling factor is calculated separately for each set of parameter combinations, where minimum popn size is 1000
#Version with reduced priors - 5k-50k.
import random
import sys
import math

#defined constants:
num_gen_change = 200
startingID = 1
endingID = 5000
N_range_lower = 5000
N_range_upper = 50000

def get_growth_factor(N_anc, N_cur):
    if N_anc > N_cur:
        N_fold = N_anc/float(N_cur)
        g = math.log(float(N_fold))/float(num_gen_change)
        growth_factor = 1.0 - g
    else:
        N_fold = N_cur/float(N_anc)
        g = math.log(float(N_fold))/float(num_gen_change)
        growth_factor = 1.0 + g
    return growth_factor

def get_scaling_factor(N_anc, N_cur):
    N_min = min(N_anc, N_cur)
    scaling_factor = float(N_min)/5000.0
    if scaling_factor > 200.0:
        scaling_factor = 200.0
    return scaling_factor

def sample_loguniform(lower_bound, upper_bound):
    Na_abs = round(math.exp(random.uniform(math.log(lower_bound), math.log(upper_bound))))
    Nc_abs = round(math.exp(random.uniform(math.log(lower_bound), math.log(upper_bound))))
    return (Na_abs, Nc_abs)

def sample_uniform(N_range_lower, N_range_upper):
    Na_abs = round(random.uniform(N_range_lower, N_range_upper))
    Nc_abs = round(random.uniform(N_range_lower, N_range_upper))
    return (Na_abs, Nc_abs)

def sample_DFE():
    f0_0 = random.randint(0, 100)
    f1_0 = random.randint(0, 100)
    f2_0 = random.randint(0, 100)
    f3_0 = random.randint(0, 100)
    f0 = round(f0_0/float(f0_0 + f1_0 + f2_0 + f3_0))
    f1 = round(f1_0/float(f0_0 + f1_0 + f2_0 + f3_0))
    f2 = round(f2_0/float(f0_0 + f1_0 + f2_0 + f3_0))
    f3 = round(f3_0/float(f0_0 + f1_0 + f2_0 + f3_0))
    return (f0, f1, f2, f3)

#Read discrete DFE from file:
f_dfe = open("/home/pjohri1/eqm_disc_parameters+5.txt", 'r')
l_dfes = []
for line in f_dfe:
    line1 = line.strip('\n')
    if "sim" in line:
        line2 = line1.split('\t')
        if line2[0] != "simID":
            s_dfe = line2[1] + "," + line2[2] + "," + line2[3] + "," + line2[4]
            l_dfes.append(s_dfe)
f_dfe.close()

result = open("/home/pjohri1/BgsDfeDemo_Human/demo_dfe_SingExon_human_logunif_parameters_v2.txt", 'w+')
result.write("simID" + '\t' + "f0" + '\t' + "f1" + '\t' + "f2" + '\t' + "f3" + '\t' +"Na" + '\t' + "Nc" + '\t' + "Na_scaled" + '\t' + "Nc_scaled" + '\t' + "scaling_factor" + '\t' + "model" + '\t' + "Nfold" + '\t' + "g_factor" + '\n')
simID = int(startingID)
dfe_num = 0
while simID <= int(endingID):
    t_N = sample_loguniform(N_range_lower, N_range_upper)
    #t_dfe = sample_DFE()
    t_dfe = l_dfes[dfe_num].split(",")
    scaling_factor = get_scaling_factor(t_N[0], t_N[1])
    Na = round(t_N[0]/float(scaling_factor))
    Nc = round(t_N[1]/float(scaling_factor))
    g_factor = get_growth_factor(Na, Nc)
    if Na > Nc:
        model = "decline"
        Nfold = "-" + str(round(Na/float(Nc),2))
    elif Na < Nc:
        model = "growth"
        Nfold = str(round(Nc/float(Na),2))
    elif Na==Nc:
        model = "eqm"
        Nfold = str(0.00)
    result.write("sim" + str(simID) + '\t' + str(t_dfe[0]) + '\t' + str(t_dfe[1]) + '\t' + str(t_dfe[2]) + '\t' + str(t_dfe[3]))
    result.write('\t' + str(t_N[0]) + '\t' + str(t_N[1]) + '\t' + str(Na) + '\t' + str(Nc) + '\t' + str(scaling_factor) + '\t' + model + '\t' + Nfold + '\t' + str(g_factor) + '\n')
    simID = simID + 1
    dfe_num += 1
    if dfe_num == len(l_dfes):
        dfe_num = 0

result.close()
print ("Finished")

