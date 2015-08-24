"""
main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
for every run, just for every project.

To create new project files save as changing description, then search < > filling in with the appropriate info.
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

#-------------------------------------------------------------------------------------------------------------------------------------
# Parameters that must be set - Paths, dates etc. These need to be setup. As more dates are added simply add these to the list here
#-------------------------------------------------------------------------------------------------------------------------------------
project_name = 'H7N9';
pathogen_name = 'Influenza'

save_path = os.path.join(ROOT_PATH, 'ArrayData', pathogen_name, project_name)
load_path = save_path


# Dates of experiments, each one is a directory where GPR files are saved.
# multiple dates can be entered, but must have a corresponding prefix
experiment_dates = ['06_03_2015'] #['08_21_2014', '08_22_2014', '08_25_2014', '09_04_2014', '09_05_2014']; 

# all slides and slideToSampleMapping file for each date have the exact same filePrefix
exp_prefix = [None]*len(experiment_dates) # initalize empty list of size experiment_dates
for i, exp in enumerate(experiment_dates):
  exp_str = exp.replace('20', '') # strip away the 20 from the year part of the date
  exp_prefix[i] = "".join([project_name, '_', exp_str])

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


#--------------------------------------------------------------------------------#
# Read in all mat files of array data, and strip them out from the matstructs 
# into a dataframe:
#-------------------------------------------------------------------------------#
arr_df = None
for i, exp in enumerate(experiment_dates):

    mapping_filename = os.path.join(load_path, exp, "".join([exp_prefix[i], 'SlideToSampleMappings.txt'])) 
    array_data_filename = os.path.join(save_path, exp, "".join([exp_prefix[i], '_arrayData_', type_flag, '.mat']))
    print(array_data_filename)
    
    d = io.loadmat(array_data_filename, struct_as_record=False)['arrayData']
    matstruct = d[0][0]
    ptids = [matstruct.ptids[0][i][0] for i in np.arange(matstruct.ptids[0].shape[0])]
    group_names = [matstruct.groupNames[0][i][0] for i in np.arange(matstruct.ptids[0].shape[0])]
    
    res = matstruct.responseMatrix[0][0]
    antigens = [matstruct.antigenNames[i][0][0] for i in np.arange(matstruct.antigenNames.shape[0])]

    columns = antigens + ['group']
    data = np.column_stack((res, np.asarray(group_names)))

    if arr_df is None:
        arr_df = pd.DataFrame(data, index=ptids, columns=columns)
    else:
        arr_df = pd.concat((arr_df, pd.DataFrame(data, index=ptids, columns=columns)), axis=0)


#-----------------------------------------------------------------#
# Background subtraction

# remove negative control responses as measured by BSA alone:
# background subtraction: subtract maximal BSA response for each antigen
#-----------------------------------------------------------------#

bsa_inds = arr_df.index.to_series().str.contains('BSA')
bg_df = arr_df[bsa_inds]
if len(bg_df.shape) == 1:
    bg_responses = np.asarray(bg_df[antigens])
else:
    bg_responses = np.asarray(bg_df[antigens].max())
fg_inds = ~bsa_inds

bg_df = arr_df[bsa_inds]
arr_df = arr_df[fg_inds]

curr_data = np.array(arr_df.as_matrix(columns=[antigens]), dtype=float) - bg_responses.T
curr_data[np.where(curr_data < 0)] = 0
arr_df[antigens] = curr_data
arr_df[antigens] = arr_df[antigens].astype(float)

#-----------------------------------------------------------------#
# Index and group dictionary setup:

# create dictionary with index sets for all HA and NA proteins
# These are antigen mask sets...
# Note that these can only be defined After creating the 
# arr_df DataFrame
#-----------------------------------------------------------------#
ind_dict = {}  # indices of cols for each strain:
for p, s in zip(prot_names, prot_strs):
    # define indices into columns of all specific HA and NA seqs:
    ind_dict[p] = [col for col in arr_df.columns if (col[:6] == p and col[7:].isdigit())]  


# create dictionay of all exp groups in the data where data is indices into the arr_df
group_inds = {}
for e in exp_groups:
    group_inds[e] = arr_df.loc[arr_df['group'] == e].index

# label dictionary from group to a numeric label
group_labels = {}
for l, e in enumerate(exp_groups):
    group_labels[e] = l

# dicionary for time-points, i.e. for pre/post etc.
time_dict = {}
# if pre/post data:
#time_dict['Pre'] = arr_df.group.str.contains('pre')
#time_dict['Post'] = arr_df.group.str.contains('post')

time_dict['all'] = [True]*arr_df.shape[0] # point to all data if not timepoints


# add numeric label column for each group
arr_df['group_label'] = arr_df['group'].map(group_labels)

#-----------------------------------------------------------------#
# Breadth and Magnitude summary statistics for array data:
#-----------------------------------------------------------------#
# compute summary statistics of the array - breadth and magnitude and store in new columns in arr_df
# uses the ind_dict from above
for p, s in zip(prot_names, prot_strs):
    arr_df[s + '_mag'] = arr_df[ind_dict[p]].sum(axis=1)  # insert new columns into dataframe for overall magnitude for each strain

    # breadth - response is positive if it is above the mean response of that antigen across all samples
    for col in ind_dict[p]:
        arr_df[col + '_binarized'] = arr_df[col].map(lambda s: 1 if s > 2000 else 0)

    arr_df[s + '_breadth'] = arr_df[[col for col in arr_df.columns if (col[:6] == p and col.endswith('_binarized'))]].sum(axis=1)

#-------------------------------------------------------------------------------#
# analysis 
#-------------------------------------------------------------------------------#
  
#-------------------------------------------------------------------------------#
# 1. Clustering analysis:
#-------------------------------------------------------------------------------#
# Cluster using Andrew's package 
# (complete linkage using Spearman correlation coefficient):

dMat = {}  # distance matrices
Z_struct = {}  # clustering struct
dend = {}  # dendrogram struct
clusters = {}  # cluster labels
cluster_treatment_stats = {}
pred_treatment_labels = {}  # predicted treatment labels based on clustering

# cluster a given timepoint (can be Post, Pre, all etc.)
post_inds = time_dict['all'] 
p_labels = np.unique(arr_df[post_inds].group_label.values)

# clustering occurs for each protein separately: i.e. H7, N9, H1, etc.
for k in ind_dict.keys():
    # use Andrew's package which allows clustering using Spearman distances (sch.linkage, and pdist do not support this for some reason, unlike Matlab)
    (dMat[k], Z_struct[k], dend[k]) = hcp.computeHCluster(arr_df[post_inds][ind_dict[k]], method='complete', metric='spearman')
    clusters[k] = sch.fcluster(Z_struct[k], t=num_clusters, criterion='maxclust')


#-------------------------------------------------------------------------------#
# 2. Statistical comparisons between groups
#-------------------------------------------------------------------------------#
  # ranksum tests 
  #[expGroupStruct(i).totalMagP] = ranksum(totalMagnitude(WT_postInds),totalMagnitude(Ob_postInds));
  #[expGroupStruct(i).H7N9magP]  = ranksum(H7N9magnitude(WT_postInds),H7N9magnitude(Ob_postInds));


# compute median responses by treatment group:
group_medians = {}
group_stds = {}
for p, s in zip(prot_names, prot_strs):
    for g in group_inds.keys():
        curr_df = arr_df.loc[group_inds[g]][ind_dict[p]]
        group_medians[g, s] = curr_df.apply(np.median, axis=0)
        group_stds[g, s] = curr_df.apply(np.std, axis=1)




