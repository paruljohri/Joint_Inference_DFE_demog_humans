README

1. Bash script to filter .ms files and calculate stats:
run_SingExon_stats_array_osg_filtered.sh

This bash script uses:
-filter_ms_phase3_data.py //to filter simulated data
-statistics_bigwindow_pylibseq_SingExon_osg_human_filtered.py //to calculate stats using pylibseq for each exon
-get_final_statistics_SingExon_filtered.R //summarize the stats across exons

2. Sort exons as high and low divergence:
sort_exons_by_div.py

3.Bash script to filter .ms files and calculate stats by different divergence categories:
run_SingExon_stats_array_osg_filtered_by_category.sh

This script uses:
-filter_ms_phase3_data.py //to filter simulated data
-statistics_bigwindow_pylibseq_SingExon_osg_human_filtered_by_category.py //to calculate stats using pylibseq for each exon
-get_final_statistics_SingExon_filtered.R //summarize the stats across exons

4. Make the final table of stats that can be used for ABC inference
make_table_statistics_SingExon_osg.py

This script uses the parameter file:
demo_dfe_SingExon_human_logunif_parameters_v2.txt //can be found in the folder - Joint_Inference_DFE_demog_humans/ABC/Step1-simulations_on_OSG/
