""" utiltlity functions for antigen array data Hertz Lab"""
import pandas as pd
import matplotlib.mlab as mlab
import scipy.stats
import numpy as np
import itertools
import os as os
from scipy import io
import hclusterplot as hcp
import scipy.cluster.hierarchy as sch
import ipdb as ipdb

__all__ = ['']

# def compare_response_peaks(resp_df, peak_threshold=30000):
#     """
#     compare the responses of peaks in the resp_df pandas dataframe (antigen array data).

#     Parameters:
#     ----------
#     resp_df: pandas.DataFrame
#         dataframe of responses of the samples to compare
#     peak_threshold: int
#         threshold defining response peaks. Default set to 30,000. 
#     Returns:
#     -------
#         output of responses of all samples to peaks in each sample in the dataframe
#     """
#     resp_df = 

def load_array_data(load_path, project_name, exp_dates, type_flag='median'):

    """
    Load array data from a given load_path, from matlab mat files for the given 
    experimental dates. Currently assumes all exp_dates are from the same 
    project (i.e. load_path). Assumes each date has the standard SlideToSampleMappings.txt file.

    Parameters:
    ----------
    load_path: string
        path of directory in which the exp_date directories are to be found.
    project_name: string
        name of current project.
    type_flag: [string | 'median' (default)]
        type of mat file: can be 'median' or 'mean'
    exp_dates: List
        list of dates from which to load data

    Returns:
    -------
    arr_df: pandas.DataFrame
        dataframe of array data
    antigens: list
        list of antigens that are included in the arr_df
    """
    # assumes all slides and slideToSampleMapping file for each date have the exact same filePrefix
    exp_prefix = [None]*len(exp_dates) # initalize empty list of size experiment_dates
    for i, exp in enumerate(exp_dates):
        exp_str = exp.replace('20', '') # strip away the 20 from the year part of the date
        exp_prefix[i] = "".join([project_name, '_', exp_str])

    arr_df = None
    for i, exp in enumerate(exp_dates):

        mapping_filename = os.path.join(load_path, exp, "".join([exp_prefix[i], 'SlideToSampleMappings.txt'])) 
        array_data_filename = os.path.join(load_path, exp, "".join([exp_prefix[i], '_arrayData_', type_flag, '.mat']))
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

    return arr_df, antigens


def backgound_subtract_array_data(arr_df, antigens, bg_name='BSA', method='max'):
    """
    remove negative control responses as measured by BSA alone.
    background subtraction: subtract maximal BSA response for each antigen.
    Data can include one or more bg_name antigens and if many are provided
    will compute thier summary statistic as defined by method.

    Parameters:
    ----------
    arr_df: pandas.DataFrame
        dataframe with array data, inlcuding bg responses.
    antigens: list
        list of antigen columns in df which need to be bg subtracted.
    bg_name: [string | 'BSA' (default)]
        name of bg_antigen to use for subtraction.
    method: ['mean', 'max']
        summary statistic used to define bg response if more than one repeat of
        the negative control is provided.

    Returns:
    -------
    arr_df: pandas.DataFrame
        modified arr_df which now does not include the bg responses and
        in which the fg responses have all been subtracted.
    bg_df: pandas.DataFrame
        dataframe of the bg_responses un-modified.
    """
    arr_df.is_copy = False  # turn of automatic warnings on setting data into original arr_df
    bsa_inds = arr_df.index.to_series().str.contains('BSA')
    bg_df = arr_df[bsa_inds]
    if len(bg_df.shape) == 1:
        bg_responses = np.asarray(bg_df[antigens])
    else:
        if method == 'max':
            bg_responses = np.asarray(bg_df[antigens].max())
        elif method == 'mean':
            bg_responses = np.asarray(bg_df[antigens].mean())
    fg_inds = ~bsa_inds

    bg_df = arr_df[bsa_inds]
    arr_df = arr_df[fg_inds]
    arr_df.is_copy = False
    curr_data = np.array(arr_df.as_matrix(columns=[antigens]), dtype=float) - bg_responses.T
    curr_data[np.where(curr_data < 0)] = 0
    arr_df[antigens] = curr_data
    arr_df[antigens] = arr_df[antigens].astype(float)

    return arr_df, bg_df


def baseline_adjust_array_data_by_subject(arr_df, time_dict, antigens,
                                          baseline_key, method='subtraction'):
    """
    Adjust for baseline responses by subject. This can be done either by subrtacting
    the baseline response, or by dividing by it.
    Parameters:
    ----------
    arr_df: pandas.DataFrame
        dataframe with array data, inlcuding bg responses.
    time_dict: dictionary
        dictionary that identifies all baseline and other timepoint responses.
    antigens: list
        list of antigen columns in df which need to be baseline adjusted.
    baseline_key: string
        key of baseline responses in the time_dict dictionary. Assumes that
        this key may also be used in the names of samples - e.g. N38_D0, or
        D0_N38.
    method: ['subtraction', 'division']
        Adjustment type. Currently two types are supported:
        'subtraction' - in which baseline responses are subtracted from all other
        timepoints and negative values are zeroed out
        'division' - compute fold increase over baseline for all additional timepoints.

    Returns:
    -------
    arr_df: pandas.DataFrame
        modified arr_df in which all non-baseline responses were adjusted for baseline.
    """
    arr_df.is_copy = False  # turn of automatic warnings on setting data into original arr_df

    for s in time_dict[baseline_key]:
        
        curr_mask = arr_df.index.to_series().str.startswith(s.replace(baseline_key, ''))
        curr_mask.loc[s] = False

        if method is 'subtraction':
            arr_df.loc[curr_mask, antigens] = arr_df[curr_mask][antigens].sub(arr_df.loc[s][antigens])
        elif method is 'division':
            arr_df.loc[curr_mask, antigens] = arr_df[curr_mask][antigens].div(arr_df.loc[s][antigens])

    return arr_df

# def baseline_adjust_data(subject_arr_df, baseline_name, antigens, method='subtraction'):
#     subject_arr_df.is_copy = False
#     for t in subject_arr_df.index:
#         if t == baseline_name:
#             continue
#         else:
#             if method is 'subtraction':
#                 subject_arr_df.loc[t][antigens] = subject_arr_df.loc[t][antigens] - \
#                                             subject_arr_df.loc[baseline_name][antigens]
#             elif method is 'division':
#                 subject_arr_df.loc[t][antigens] = subject_arr_df.loc[t][antigens] / \
#                                             subject_arr_df.loc[baseline_name][antigens]

#     return subject_arr_df


def initialize_indexing_dictionaries(arr_df, prot_names, prot_strs, exp_groups, pre_post_data=False):
    """
    Index and group dictionary setup:
    create dictionary with index sets for all HA and NA proteins
    These are antigen mask sets... 
    
    Parameters:
    ----------
    arr_df: pandas.DataFrame
        dataframe of array data
    prot_names: List
        list of strings of protein antigens from which peptides are on the array_data_filename
        defined in init_array_data_params_<project>
    prot_strs: List
        names of proteins in prot_names used for plotting and vizualization.
        defined in init_array_data_params_<project>
    exp_groups: list
        names of experimental groups
        defined in init_array_data_params_<project>
    pre_post_data: bool 
        if True data contains longitudinal pre/post timepoints. default set to false.
        Used for time_dict dictionary
    
    Returns:
        ind_dict: dictionary
            dictionary mapping prot_names to columns in the arr_df
        group_inds: dictionary
            dictionary mapping group names to rows in the arr_df dataframe
        time_dict: dictionary
            dictionary separating samples by timepoints (defaults to all data).
            Note that always includes an 'all' key which includes all data.
        group_labels: dictionary
            dictionary mapping group names to numerical group labels
    """
    ind_dict = {}  # indices of cols for each strain:
    for p, s in zip(prot_names, prot_strs):
        # define indices into columns of all specific HA and NA seqs:
        ind_dict[p] = [col for col in arr_df.columns if (col[:6] == p and col[7:].isdigit())]  


    # create dictionay of all exp groups in the data where data is indices into the arr_df
    group_inds = {}
    for e in exp_groups:
        group_inds[e] = arr_df.loc[arr_df['group'] == e].index


    # dicionary for time-points, i.e. for pre/post etc.
    time_dict = {}
    if pre_post_data:
        time_dict['Pre'] = arr_df.group.str.contains('pre')
        time_dict['Post'] = arr_df.group.str.contains('post')

    # all data is always included in this dictionary for default.
    time_dict['all'] = [True]*arr_df.shape[0] # point to all data if not timepoints

     # label dictionary from group to a numeric label
    group_labels = {}
    for l, e in enumerate(exp_groups):
        group_labels[e] = l

    return ind_dict, time_dict, group_inds, group_labels


def add_numeric_labels_to_arr_df(arr_df, group_labels):
    """add numeric label column for each group"""
    arr_df['group_label'] = arr_df['group'].map(group_labels)
    return arr_df


def compute_breadth_and_magnitude_summary_stats(arr_df, prot_names, prot_strs, ind_dict):
    """
    Breadth and Magnitude summary statistics for array data:
    compute summary statistics of the array - breadth and magnitude and stores them in 
    new columns in arr_df. Uses the ind_dict dictionary for each protein in prot_names.
    
    Parameters:
    ----------
    arr_df: pandas.DataFrame
        dataframe of array data
    prot_names: List
        list of strings of protein antigens from which peptides are on the array_data_filename
    prot_strs: List
        names of proteins in prot_names used for plotting and vizualization.
    ind_dict: dictionary
        dictionary mapping prot_names to columns in the arr_df    

    Returns:
    -------
    arr_df: pandas.DataFrame
        dataframe with magnitude and breadth columns added for each protein in prot_strs.
    """
    for p, s in zip(prot_names, prot_strs):
        # insert new columns into dataframe for overall magnitude for each strain
        arr_df[s + '_mag'] = arr_df[ind_dict[p]].sum(axis=1)  

    # breadth - response is positive if it is above the mean response of that antigen across all samples
    for col in ind_dict[p]:
        arr_df[col + '_binarized'] = arr_df[col].map(lambda s: 1 if s > 2000 else 0)

    arr_df[s + '_breadth'] = \
        arr_df[[col for col in arr_df.columns if (col[:6] == p and col.endswith('_binarized'))]].sum(axis=1)

    return arr_df

def cluster_array_data_by_proteins(arr_df, sample_inds, num_clusters, ind_dict, method='complete', metric='spearman'):
    """
    Cluster array data for each protein antigen separately using linkage clustering

    Parameters:
    ----------
    arr_df: pandas.DataFrame
        dataframe of array data
    sample_inds: list
        list of samples that we want to cluster - allows clustering a subset of the data in arr_df.
    num_clusters: int
        number of clusters - arbitrary parmeter set by the user.
    ind_dict: dictionary
        dictionary mapping prot_names to columns in the arr_df
    method: string
        linkage clustering method. Can be 'average', 'single' and 'complete' (default).
    metric: string
        distance metric used for comparing response vectors. Default set to 'spearman' (rank order correlation).

    Returns:
    -------
    dist_mat: dictionary
        dictionary of distance matrices whose keys are the keys in ind_dict (proteins)
    Z_struct: dictionary
        clustering structure (matlab style) returned by linkage indexed by ind_dict.keys()
    dend: dictionary
        clustering dendrograms for each index set in ind_dict
    clusters: dictionary
        cluster assignment of each datapoint indexed by ind_dict.keys()
    """
    dist_mat = {}  # distance matrices
    Z_struct = {}  # clustering struct
    dend = {}  # dendrogram struct
    clusters = {}  # cluster labels
    
    # clustering occurs for each protein separately: i.e. H7, N9, H1, etc.
    for k in ind_dict.keys():
        # use Andrew's package which allows clustering using Spearman distances 
        # (sch.linkage, and pdist do not support this for some reason, unlike Matlab)
        (dist_mat[k], Z_struct[k], dend[k]) = \
            hcp.computeHCluster(arr_df[sample_inds][ind_dict[k]], method='complete', metric='spearman')
        clusters[k] = sch.fcluster(Z_struct[k], t=num_clusters, criterion='maxclust')

    return dist_mat, Z_struct, dend, clusters


def compute_exp_group_summary_stats_by_protein(arr_df, prot_names, prot_strs, group_inds, ind_dict):
    """ 
    compute median and std of the responses of each expereimental group (defined by group_inds)
    to each ot the proteins in prot_names.

    Parameters:
    ----------
    arr_df: pandas.DataFrame
        dataframe of array data
    prot_names: list
        list of strings of protein antigens from which peptides are on the array_data_filename
    prot_strs: list
        names of proteins in prot_names used for plotting and vizualization.
    group_inds: dictionary
        dictionary mapping group names to rows in the arr_df dataframe
    ind_dict: dictionary
        dictionary mapping prot_names to columns in the arr_df

    Returns:
    -------
    group_medians: dictionary
        median response of each group to each individual antigen
        for a given protein in prot_strs
    group_stds: dictionary
        standard deviation of the responses of each group to each individual antigen
        for a given protein in prot_strs
    """
    group_medians = {}
    group_stds = {}
    for p, s in zip(prot_names, prot_strs):
        for g in group_inds.keys():
            curr_df = arr_df.loc[group_inds[g]][ind_dict[p]]
            group_medians[g, s] = curr_df.apply(np.median, axis=0)
            group_stds[g, s] = curr_df.apply(np.std, axis=1)

    return group_medians, group_stds

