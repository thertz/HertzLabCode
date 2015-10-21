""" Template script for analyzing array data from a specific project"""
from __future__ import division, print_function
import sys as sys
sys.path.append('../Utils')
import os as os
import numpy as np
import scipy.stats
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
project_name = 'H7N9';

# Names of HA and NA proteins that are on the array. 
prot_names = ['SHA_ha', 'SHA_na', 'Cal_ha', 'Cal_na']
prot_strs = ['H7', 'N9', 'H1', 'N1']  # The names for figure plotting

# Labels of the experimental groups (not including BSA, or any other background):
exp_groups = ['Normal', 'Obese']
exp_group_prefixes = ['Obese_', 'Normal_']

# Summary statistics (breadth and magnitude) that we would like to visualize. 
# if set to empty, will plot all.
arr_summary_stats = ['H7_mag', 'N9_mag', 'H1_mag', 'N1_mag']

# Paths:
ROOT_PATH = '/Users/thertz/Dropbox/HertzLab/'
load_path = os.path.join(ROOT_PATH, 'ArrayData', pathogen_name, project_name)
fig_path = load_path + '/TestFigs/'


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
experiment_dates = ['06_03_2015']
arr_df, antigens = amutils.load_array_data(load_path, project_name, experiment_dates)


"""Background subtraction:"""
arr_df, bg_df = amutils.backgound_subtract_array_data(arr_df, antigens, bg_name='BSA')


"""Index and group dictionary setup:"""
ind_dict, time_dict, group_inds, group_labels = \
    amutils.initialize_indexing_dictionaries(arr_df, prot_names, prot_strs, exp_groups)

"""add numeric label column to arr_df for stat testing."""
arr_df = amutils.add_numeric_labels_to_arr_df(arr_df, group_labels)

"""Breadth and Magnitude summary statistics of the array data:"""
arr_df = amutils.compute_breadth_and_magnitude_summary_stats(arr_df, prot_names, prot_strs, ind_dict)


"""
Clustering analysis:

Cluster using Andrew's package (complete linkage using Spearman correlation coefficient):
"""

"""Unsupervised clustering parameters:"""
num_clusters = 4 # number of clusters in the data (or expected number)
# minimal threshold for responces to be included in clusterSamplesByResponseVectors, point below data is noise
minResponseThreshold = 2000 

dist_mat, Z_struct, dend, clusters = \
    amutils.cluster_array_data_by_proteins(arr_df=arr_df, sample_inds=time_dict['all'],
                                           num_clusters=num_clusters, ind_dict=ind_dict)


"""Statistical comparisons between groups"""
group_medians, group_stds = \
    amutils.compute_exp_group_summary_stats_by_protein(arr_df, prot_names, prot_strs, group_inds, ind_dict)


"""Data Visualzation:"""


"""Responses by groups:"""
"""
1. Median responses by group:
Plot figures of median responses by treatment group - Note here names are from prot_strs and NOT prot_names
"""
for p in ['H7',  'N9']:
    amp.plot_median_responses_by_exp_groups(group_medians, exp_groups, prot_str=p,
                                            fig_path=None, fig_prefix=None, fig_size=(18,11), y_lims=[0, 30000])
    plt.show()
    display(HTML("<HR>"))


"""
2. Raw responses by group:
Plot responses by groups for each protein separately - groups can be prot_strs, or any subset of these:
the second list is prot_strs for figure titles. 
you can use this for all proteins:
for p, s in zip(prot_names, prot_strs):
"""
for p, s in zip(['SHA_ha', 'SHA_na'], ['H7', 'N9']):
    amp.plot_responses_by_exp_groups(arr_df=arr_df, antigen_inds = ind_dict[p],
                                     exp_groups=exp_groups, fig_path=None, fig_prefix=s, y_lims=y_lims)

"""
3. Boxplot of summary stats by group:
plot boxplots of all groups - magnitude and breadth for each protein 
"""
amp.plot_summary_stat_boxplots_by_exp_groups(arr_df, arr_summary_stats, sample_inds=None, fig_path=None)
    


"""Clustering plots"""

"""
1. Dendrograms:
"""
amp.plot_clustering_dendrograms(Z_struct=Z_struct, prot_names=prot_names, labels=arr_df.index, fig_path=None)


"""
2. Raw responses by cluster:
"""
for p, s in zip(['SHA_ha', 'SHA_na'], ['H7', 'N9']):
    amp.plot_raw_responses_by_clusters(arr_df, ind_dict[p], num_clusters, clusters=clusters[p],
                                          fig_path=None, fig_prefix=s, fig_size=(18,11), y_lims=[0, 30000])
 
"""
3. Median responses by cluster:
"""
for p, s in zip(['SHA_ha', 'SHA_na'], ['H7', 'N9']):
    amp.plot_median_responses_by_clusters(arr_df, ind_dict[p], num_clusters, clusters=clusters[p],
                                          fig_path=None, fig_prefix=s, fig_size=(18,11), y_lims=[0, 30000])

"""
4. Summary stat boxplots by clusters:
"""
amp.plot_summary_stat_boxplots_by_clusters(arr_df, clusters, prot_names, arr_summary_stats,
                                           sample_inds=None, fig_path=None)