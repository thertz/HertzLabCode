""" Template script for analyzing array data from a specific project"""
from __future__ import division, print_function
import sys as sys
sys.path.append('../Utils')
import os as os
import numpy as np
import pandas as pd
import scipy.stats
import itertools
import matplotlib.pyplot as plt
import hclusterplot as hcp
import myboxplot as mbp
import scipy.cluster.hierarchy as sch
import statsmodels.api as sm
from IPython.display import display, HTML
import amplotlib as amp
import amutils as amutils

"""Project specific parameters that must be set:"""
# These are all just examples that can and should be modified.
pathogen_name = 'Influenza'
project_name = 'IDRI';

# Names of HA and NA proteins that are on the array. 
prot_names = ['VN1203', 'VN1194', 'Thai16', 'Kan-1', 'r0405050', 'HK483', 'Indo05',
              'Anhui', 'Turkey', 'Monoglia', 'Duck-Hunan']
prot_strs = prot_names  # The names for figure plotting

# Labels of the experimental groups (not including BSA, or any other background):
exp_groups = ['D0', 'D21', 'D42']
#exp_group_prefixes = ['Obese_', 'Normal_'] # no prefixes here.
exp_group_suffixes = ['_D0', '_D21', '_D42']

# Summary statistics (breadth and magnitude) that we would like to visualize. 
# if set to empty, will plot all.
arr_summary_stats = ['VN1203_mag', 'Indo05_mag']

# Paths:
ROOT_PATH = '/Users/thertz/Dropbox/HertzLab/'
# ROOT_PATH = 'c:/Users/tomer.hertz/Dropbox/HertzLab/'

load_path = os.path.join(ROOT_PATH, 'ArrayData', pathogen_name, project_name)
fig_path = load_path + '/Figs/'


"""Additional parameters that have usual defaults:"""
type_flag  = 'median' # uses the median over replicates of an antigen. Can also be 'mean'
color_flags = ['635']
color_tags = ['IgG']
font_size = 12

y_lims = [0, 65000] # max height of y in graphs
y_lims_summary = [0, 32500] # max height of y in summary stat graphs

# Threshold by which to flag antigens as ones that have high background and should be removed from analysis 
# (typically using a BSA negative control).
bg_threshold = 2000 

# summary stat used for group comparisons and for plotting gropu responses - and cluster responses
summary_method = 'median' 


"""
    Read in data from matlab mat files
    This is obsolete and will be modified in future version.
"""
experiment_dates = ['05_08_2015', '05_13_2015', '05_15_2015']
arr_df, antigens = amutils.load_array_data(load_path, project_name, experiment_dates)

"""Background subtraction:"""
arr_df, bg_df = amutils.backgound_subtract_array_data(arr_df, antigens, bg_name='BSA')

# small tweak to unify timepoint labels (D0 with d0)
arr_df.index = arr_df.index.to_series().str.upper()

"""
Index and group dictionary setup:
Here we had multiple strains made with overlapping peptides, so that some peptides match several strains.
Therefore we need to create dictionaries that allow overalp of some peptides to several strains.
"""
# parse IDRI peptide set cheat sheet ot get correct naming of peptides
idri_pep_df = pd.read_excel(load_path + '/IDRI_peptide_mapping.xlsx')

# index dictionaries - here some peptides belong to multiple strains
ind_dict = {}
beg_ind_dict = {}
for p in prot_names:
    curr_df = idri_pep_df[idri_pep_df[p] == 1]
    ind_dict[p] = curr_df.pepName
    beg_ind_dict[p] = curr_df.Position

# create dictionay of all exp groups in the data where data is indices into the arr_df
group_inds = {}
for e in exp_groups:
    group_inds[e] = arr_df.loc[arr_df['group'] == e].index

# label dictionary from group to a numeric label
group_labels = {}
for l, e in enumerate(exp_groups):
    group_labels[e] = l

# currently groups are timepoints
time_dict = {} # currently exp_groups are timepoints!
for e in exp_groups:
    time_dict[e] = arr_df['group'] == e
time_dict['all'] = [True]*arr_df.shape[0] # point to all data if not timepoints 


"""add numeric label column to arr_df for stat testing."""
arr_df = amutils.add_numeric_labels_to_arr_df(arr_df, group_labels)

"""Baseline adjustment - D0 responses here"""
# first turn all response values to antigens into floats!
all_antigens = np.unique([item for sublist in ind_dict.values() for item in sublist])
arr_df[antigens] = arr_df[antigens].astype('float')

# now adjust for baseline response of each sample:
raw_arr_df = arr_df.copy()

# note that here gropu_inds is used at the time frame... (NEEDS TO BE THOUGHT ABOUT)
arr_df = amutils.baseline_adjust_array_data_by_subject(arr_df, time_dict=group_inds, antigens=all_antigens,
                                                       baseline_key='D0', method='subtraction')
# zero out all negative rsponses
arr_df[antigens] = arr_df[antigens].clip(0, None)

"""Breadth and Magnitude summary statistics of the array data:"""
arr_df = amutils.compute_breadth_and_magnitude_summary_stats(arr_df, prot_names, prot_strs, ind_dict)
raw_arr_df = amutils.compute_breadth_and_magnitude_summary_stats(raw_arr_df, prot_names, prot_strs, ind_dict)

"""
Clustering analysis:

Cluster using Andrew's package (complete linkage using Spearman correlation coefficient):
"""

"""Unsupervised clustering parameters:"""
num_clusters = 4 # number of clusters in the data (or expected number)
# minimal threshold for responces to be included in clusterSamplesByResponseVectors, point below data is noise
minResponseThreshold = 1000 

# dist_mat, Z_struct, dend, clusters = \
#     amutils.cluster_array_data_by_proteins(arr_df=arr_df, sample_inds=time_dict['D21'],
#                                            num_clusters=num_clusters, ind_dict=ind_dict)


dist_mat, Z_struct, dend, clusters = \
    amutils.cluster_array_data_by_proteins(arr_df=raw_arr_df, sample_inds=time_dict['all'],
                                           num_clusters=num_clusters, ind_dict=ind_dict)


"""Statistical comparisons between groups"""
group_medians, group_stds = \
    amutils.compute_exp_group_summary_stats_by_protein(arr_df, prot_names, prot_strs, group_inds, ind_dict)

raw_group_medians, raw_group_stds = \
    amutils.compute_exp_group_summary_stats_by_protein(raw_arr_df, prot_names, prot_strs, group_inds, ind_dict)


"""Data Visualzation:"""


"""Responses by groups:"""
"""
1. Median responses by group:
Plot figures of median responses by treatment group - Note here names are from prot_strs and NOT prot_names
"""
for p in ['VN1203',  'Indo05']:
    amp.plot_median_responses_by_exp_groups(raw_group_medians, exp_groups, prot_str=p,
                                            fig_path=fig_path, fig_prefix='Raw_', fig_size=(18,11), y_lims=[0, 50000])


for p in ['VN1203',  'Indo05']:
    amp.plot_median_responses_by_exp_groups(group_medians, ['D21', 'D42'], prot_str=p,
                                            fig_path=fig_path, fig_prefix=None, fig_size=(18,11), y_lims=[0, 10000])



"""
2. Raw responses by group:
Plot responses by groups for each protein separately - groups can be prot_strs, or any subset of these:
the second list is prot_strs for figure titles. 
you can use this for all proteins:
for p, s in zip(prot_names, prot_strs):
"""
for p in ['VN1203', 'Indo05']:
    amp.plot_responses_by_exp_groups(arr_df=raw_arr_df, antigen_inds = ind_dict[p],
                                     exp_groups=exp_groups, fig_path=fig_path, fig_prefix=p + '_Raw', y_lims=y_lims)

for p in ['VN1203', 'Indo05']:
    amp.plot_responses_by_exp_groups(arr_df=arr_df, antigen_inds = ind_dict[p],
                                     exp_groups=['D21', 'D42'], fig_path=fig_path, fig_prefix=p, y_lims=y_lims)


"""
3. Boxplot of summary stats by group:
plot boxplots of all groups - magnitude and breadth for each protein 
"""
amp.plot_summary_stat_boxplots_by_exp_groups(raw_arr_df, arr_summary_stats, sample_inds=None, fig_path=fig_path)
    

"""Clustering plots"""

"""
1. Dendrograms:
"""
#amp.plot_clustering_dendrograms(Z_struct=Z_struct, prot_names=prot_names, labels=arr_df.index[time_dict['D21']], fig_path=None)
amp.plot_clustering_dendrograms(Z_struct=Z_struct, prot_names=['VN1203', 'Indo05'], labels=arr_df.index, fig_path=fig_path, fig_prefix='Raw_')

"""
2. Raw responses by cluster:
"""
for p, s in zip(['VN1203', 'Indo05'], ['VN1203', 'Indo05']):
    amp.plot_raw_responses_by_clusters(raw_arr_df, ind_dict[p], num_clusters, clusters=clusters[p],
                                          fig_path=fig_path, fig_prefix=s + '_Raw', fig_size=(18,11), y_lims=[0, 50000])
 
"""
3. Median responses by cluster:
"""
for p, s in zip(['VN1203', 'Indo05'], ['VN1203', 'Indo05']):
    amp.plot_median_responses_by_clusters(raw_arr_df, ind_dict[p], num_clusters, clusters=clusters[p],
                                          fig_path=fig_path, fig_prefix=s + '_Raw', fig_size=(18,11), y_lims=y_lims)

"""
4. Summary stat boxplots by clusters:
"""
amp.plot_summary_stat_boxplots_by_clusters(raw_arr_df, clusters, ['VN1203', 'Indo05'], arr_summary_stats,
                                           sample_inds=None, fig_path=fig_path)