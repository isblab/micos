import numpy as np
import pandas as pd
import math
import glob
import sys
import os

# Change this to the location of your PMI_analysis folder
sys.path.append('/home/muskaan/imp-clean/PMI_analysis/pyext/src')
import analysis_trajectories

#################################
########### MAIN ################
#################################

nproc = 50
top_dir = sys.argv[1]
analys_dir = './model_analysis/'

# Check if analysis dir exists
if not os.path.isdir(analys_dir):
    os.makedirs(analys_dir)

# How are the trajectories dir names
dir_head = 'run_'
out_dirs = glob.glob(top_dir+'/'+dir_head+'*/')
# print(out_dirs)
################################
# Get and organize fields for
# analysis
################################
# Read the total score, plot
# and check for score convengence
XLs_cutoffs = {'BDP_PIR_human':52.0, 'BDP_PIR_mouse':52.0, 'DSS_BS3':35.0, 'DHSO_DSSO':30, 'DSSO_mouse':30 }

# Load module
AT = analysis_trajectories.AnalysisTrajectories(out_dirs,
                          dir_name=dir_head,
                          analysis_dir = analys_dir,
                          nproc=nproc)

# Define restraints to analyze
AT.set_analyze_XLs_restraint(XLs_cutoffs = XLs_cutoffs,
				ambiguous_XLs_restraint = True)
AT.set_analyze_Connectivity_restraint()
AT.set_analyze_Excluded_volume_restraint()
# AT.set_analyze_SurfaceLocalization_restraint()
# AT.set_analyze_score_only_restraint(handle='SurfaceLocalizationRestraint', short_name='SLR', do_sum=True)
AT.set_analyze_score_only_restraint(handle='IMSLocalizationRestraint', short_name='ILR', do_sum=True)
AT.set_analyze_score_only_restraint(handle='MatrixLocalizationRestraint', short_name='MLR', do_sum=True)
AT.set_analyze_score_only_restraint(handle='ZAxialRestraint', short_name='ZAR', do_sum=True)
AT.set_analyze_score_only_restraint(handle='TransMembraneRestraint', short_name='TMR', do_sum=True)
AT.set_analyze_score_only_restraint(handle='MinimumPairDistanceBindingRestraint', short_name='MPDBR', do_sum=True)

# Read stat files
AT.read_stat_files()
AT.write_models_info()
AT.get_psi_stats()

# What scores do we cluster on?
AT.hdbscan_clustering(['EV_sum', 'XLs_sum', 'ILR_sum', 'ZAR_sum', 'TMR_sum', 'MLR_sum', 'MPDBR_sum'], min_cluster_size = 150, min_samples = 5, skip = 1) #'SLR_sum', 
AT.summarize_XLs_info(ambiguous_XLs_restraint = True)

exit()
