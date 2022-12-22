#This is to make a final table of all statistics:
#python make_table_statistics_SingExon_osg.py demo_dfe_SingExon_human_v2_stats 1000
#python make_table_statistics_SingExon_osg.py demo_dfe_SingExon_human_v2_stats_filtered_5p 2000
#python make_table_statistics_SingExon_osg.py demo_dfe_SingExon_human_v2_stats_filtered_3p 2000
#python make_table_statistics_SingExon_osg.py demo_dfe_SingExon_human_v2_stats_filtered_5p_high_div 2000
#python make_table_statistics_SingExon_osg.py demo_dfe_SingExon_human_v2_stats_filtered_5p_low_div 2000
#python make_table_statistics_SingExon_osg.py demo_dfe_SingExon_human_v2_stats_filtered_3p_high_div 2000
#python make_table_statistics_SingExon_osg.py demo_dfe_SingExon_human_v2_stats_filtered_3p_low_div 2000

import sys
folder = sys.argv[1]
tot_numsims = int(sys.argv[2])
outputfile = folder.replace("/", "_")

#Going through stats files:
result = open("/scratch/pjohri1/BgsDfeDemo_Human/" + outputfile + "_abc.txt", 'w+')

l_stats = ["thetapi_m","thetaw_m","thetah_m", "hprime_m", "tajimasd_m", "numSing_m", "hapdiv_m", "rsq_m", "D_m", "Dprime_m", "div_m", "thetapi_sd", "thetaw_sd", "thetah_sd", "hprime_sd", "tajimasd_sd", "numSing_sd", "hapdiv_sd", "rsq_sd", "D_sd", "Dprime_sd", "div_sd"]
result.write("simID" + '\t' + "f0" + '\t' + "f1" + '\t' + "f2" + '\t' + "f3" + '\t' + "Na" + '\t' + "Nc" + '\t' + "Na_scaled" + '\t' + "Nc_scaled" + '\t' + "scaling_factor" + '\t' + "model" + '\t' + "Nfold" + '\t' + "g_factor")
#Na      Nc      Na_scaled       Nc_scaled scaling_factor  model   Nfold   g_factor
for stat in l_stats:
	result.write('\t' + "func_" + stat + '\t' + "linked1_" + stat + '\t' + "linked2_" + stat)
result.write('\n')

#get parameters:
d_para = {}
f_para = open("/home/pjohri1/BgsDfeDemo_Human/demo_dfe_SingExon_human_logunif_parameters_v2.txt", 'r')
for line in f_para:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    d_para[line2[0]] = line1
f_para.close()


i = 1
mark = 0
while i <= tot_numsims:
    simID = "sim" + str(i)
    s_para = ""
	#find out if all required files exist first, for this simID:
    s_50 = ""
    #if 1 > 0:
    try:
        f = open("/scratch/pjohri1/BgsDfeDemo_Human/" + folder + "/" + simID + ".bigwinsummary", 'r')
        #read the parameter file:
        s_para = d_para[simID]
        print("para: " + s_para)

		#write numbp50 stats:
        s_50 = s_para + '\t'
        for line in f:
            line1 = line.strip('\n')
            line2 = line1.split('\t')
            if line2[0] != "func":
                for x in line2[1:]:
                    if x != "NA":
                        s_50 = s_50 + '\t' + str("{:.5f}".format(float(x)))
                    else:
                        s_50 = s_50 + '\t' + x
        f.close()
        if s_50 != "":
            result.write(s_50 + '\n')
    #else:
    except:
        print (" File not found for:" + simID)
    i = i + 1

result.close()

print ("done")



