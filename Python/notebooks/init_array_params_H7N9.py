"""
Initialize array parameters for H7N9.
"""
from __future__ import division, print_function
from scipy import io
import os as os
import pandas as pd
import numpy as np
import scipy
import scipy.stats
import hclusterplot as hcp
import myboxplot as mbp
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import statsmodels.api as sm
import amplotlib as amp
import amutils as amutils

# General parameters
cluster_plot_flag = False
plot_flag = False
save_flag = False
plt.close('all') # close all figures

#-------------------------------------------------------------------------------------------------------------------------------------
# project specific parameters - this is for H7N9, not relevant for other projects
#-------------------------------------------------------------------------------------------------------------------------------------
ROOT_PATH = '/Users/thertz/Dropbox/HertzLab/'
SAVE_PATH = ROOT_PATH + 'ArrayData/Influenza/H7N9/'
FIG_PATH = ROOT_PATH + 'ArrayData/Influenza/H7N9/AnatProject/'

# names of HA and NA proteins that are on the array. the prot_strs are the names for figure plotting
prot_names = ['SHA_ha', 'SHA_na', 'Cal_ha', 'Cal_na']
prot_strs = ['H7', 'N9', 'H1', 'N1']

arr_summary_stats = ['H7_mag', 'N9_mag', 'H1_mag', 'N1_mag']

# labels of experimental groups
#exp_groups    = ['Vac', 'AS03', 'MF59', 'PBS']; # types of adjuvants used on both Obese and WT mice
#exp_group_prefixes = ['WT_pre_', 'Ob_pre_', 'WT_post_', 'Ob_post_']
 
exp_groups = ['Normal', 'Obese']
exp_group_prefixes = ['Obese_', 'Normal_', 'BSA']

save_path = os.path.join(ROOT_PATH, 'ArrayData', pathogen_name, project_name)
load_path = save_path

#-------------------------------------------------------------------------------
# Additional parameters - project specific
#-------------------------------------------------------------------------------
num_arrays = 2;  # number of arrays on each slide.
type_flag  = 'median' # uses the median over replicates of an antigen. Can also be 'mean'
color_flags = ['635']
color_tags = ['IgG']
font_size = 12

y_lims = [0, 65000] # max height of y in graphs
y_lims_summary = [0, 32500] # max height of y in summary stat graphs
min_threshold = 10000 # minimal threshold for peak responses, used by findPeaks.
bg_threshold = 5000 # threshold by which to flag antigens as ones that have high background and should be removed from analysis (typically using a BSA negative control).
summary_method = 'median' # summary stat used for group comparisons and for plotting gropu responses - and cluster responses

# Unsupervised clustering parameters:
num_clusters = 4 # number of clusters in the data (or expected number)
#minResponseThreshold = 2000 # minimal threshold for responces to be included in clusterSamplesByResponseVectors, point below data is noise

#-------------------------------------------------------------------------------#
# Stat comparison parameters: (also project specific)
#-------------------------------------------------------------------------------#
# Filters used for treatment-blinded filtering of responses to single antigens used for statistical analysis.
# currently only one filter is implemented which is # responders within the entire treatment set (ignoring all controls) 
filterNames = ['PercentResponders']
filterThresholds = [0.2]
treatmentGroups = ['Ob_post_Vac', 'WT_post_Vac'] # Names of the two groups that are to be compared statistically below (currently supports only pairwise comparisons)

# data Labels
maskLabels = ['Cal_HA_Inds', 'Cal_NA_Inds', 'Sha_HA_Inds', 'Sha_NA_Inds']
treatmentLabels = ['WT-pre-vac', 'obese-pre-vac', 'wt-post-vac', 'obese-post-vac']


