""" Plotting library for Hertz Lab antigen array data """
import palettable
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
import matplotlib.mlab as mlab
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch
import scipy.stats
import numpy as np
import matplotlib as mpl
import seaborn as sbn
import pylab
import myboxplot as mbp
import itertools
import ipdb as ipdb

__all__ = ['plot_linear_corr', 'plot_longitudinal_responses_by_ptid', 'plot_responses_by_clusters',
           'plot_responses_by_groups', 'plot_responses_by_ptids', 'plot_2d_scatter_plot',
           'plot_reproducibility_plots']


def plot_linear_corr(x, y, x_label=None, y_label=None, title_str=None):
    """
    plot linear correlation of x vs. y, including linear regression line overlaid

    Parameters:
    ----------
    x: np.ndarray
        first data vector
    y: np.ndarray
        second data vector
    x_label: string
        label of x axis
    y_label: string
        label of y axis
    title_str: string | None (default)
        title of plot.
    
    Returns:
    -------
    f : figure handle
    """
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    line = slope * x + intercept
    f = plt.figure()
    plt.plot(x, line, 'r-', x, y, 'o')
    if x_label is not None:
        plt.xlabel(x_label)
    if y_label is not None:
        plt.ylabel(y_label)
    if title_str is None:
        plt.title("r = {:.2f}, p = {:.4f}".format(r_value, p_value))
    else:
        plt.title(title_str + " r = {:.2f}, p = {:.4f}".format(r_value, p_value))
    return f


def plot_responses_by_ptids(resp_mat, ptid_labels=None, plot_type='line', subplot_flag=False, y_lim=60000):
    """
    plot responses of a set of ptids. Can be used for reproducibility comparisons.
    
    Parameters:
    ----------
    resp_mat : np.ndarray, or pd.DataFrame
        antigen array responses where each row corresponds to responses for a single ptid
    ptid_labels : [ List | None (default)]
        ptid ids for plotting (optional)
    plot_type : ['line' (default) | 'bar']
        type of plot used to present data
    subplot_flag : [True | False]
        if True plot each ptid in a different subplot, else plot all on a single plot.

    Returns:
    -------
    f : figure handle

    """
    N, num_antigens  = resp_mat.shape
   
    width = 0.3      # the width of the bars
    num_colors = max(3,N)
    cmap_name = "palettable.colorbrewer.qualitative." + 'Set1' + "_" + str(num_colors) + ".mpl_colors"
    colors = eval(cmap_name)

    if isinstance(resp_mat, pd.DataFrame):
        resp_mat = resp_mat.as_matrix()

    if subplot_flag:
        num_plots = N
    else:
        num_plots = 1
    f, axarr = plt.subplots(num_plots, 1)
    f.set_tight_layout(True)

    curr_inds = np.arange(1, num_antigens+1)
    for i in np.arange(N):
        
        if subplot_flag:
            curr_ax = axarr[i]
        else:
            curr_ax = axarr
            if (plot_type == 'bar') and (i > 0):
                curr_inds = curr_inds + width

        if plot_type == 'line':        
            curr_ax.plot(curr_inds, resp_mat[i],color=colors[i], linestyle='-')
            
        elif plot_type == 'bar':
            curr_ax.bar(curr_inds, resp_mat[i], width=width, color=colors[i])            
    
        if (ptid_labels is not None) and subplot_flag:
            curr_ax.set_title("Response for ptid " + ptid_labels[i])

        #curr_ax.set_yticks([])
        curr_ax.set_ylim(0, y_lim)
    
    if (ptid_labels is not None) and not subplot_flag:
        plt.legend(ptid_labels)
        print(ptid_labels)
    return f


def plot_2d_scatter_plot(resp1, resp2, labels=None, title_str=None):
    """
    Plot scatter plot of responses of resp1 vs. resp2 and reorts linear
    correlation coefficient betwen the two.
 
    Parameters:
    ----------
    resp1: [np.ndarray | pd.DataFrame]
        responses of first sample
    resp2: [np.ndarray | pd.DataFrame]
        responses of second sample
    labels: [tuple | None (default)]
        labels of the two responses for axis labeling. order is (x_label, y_label)
    title_str: [string | None (default)]
        title of plot. note that title will always appear with Pearson 
        r and p reported.
    Returns:
    -------
    f : figure handle
    """
    r, p = scipy.stats.pearsonr(resp1, resp2)
    
    f, axarr = plt.subplots(1)
    f.set_tight_layout(True)
    axarr.plot(resp1, resp2, 'o')
    axarr.set_ylim(-1000, 60000)
    axarr.set_xlim(-1000, 60000)  
    if labels is not None:
        plt.xlabel(labels[0])
        plt.ylabel(labels[1])
    if title_str is None:
        plt.title("r = {:.2f}, p = {:.4f}".format(r, p))
    else:
        plt.title(title_str + " r = {:.2f}, p = {:.4f}".format(r, p))
    return f
    

def plot_reproducibility_plots(resp_mat, ptid_labels=None, title_str=None,
                               subplot_flag=False, plot_type='line'):
    """
    Plot two plots for assesing reproducibility.

    Parameters:
    ----------
    resp_mat : np.ndarray, or pd.DataFrame
        antigen array responses where each row corresponds to responses for a single ptid
        assumes only 2 rows are proivded, one for each sample.
    ptid_labels : [ List | None (default)]
        ptid ids for plotting (optional)
    title_str: [string | None (default)]
        title of plot. note that title will always appear with Pearson 
        r and p reported.
    subplot_flag : [True | False]
        if True plot response of each ptid in a different subplot, else plot all on a single plot.
    plot_type : ['line' (default) | 'bar']
        type of plot used to present data when plottting individul responses.

    Returns:
    -------
    fig_handles: list
        list of figure handles
    """
    fig_handles = []
    fig_handles.append(plot_responses_by_ptids(resp_mat=resp_mat, ptid_labels=ptid_labels,
                                               subplot_flag=subplot_flag, plot_type=plot_type))
    if isinstance(resp_mat, pd.DataFrame):
        resp1 = resp_mat.iloc[0]
        resp2 = resp_mat.iloc[1]
    else:
        resp1 = resp_mat[0]
        resp2 = resp_mat[1]

    fig_handles.append(plot_2d_scatter_plot(resp1=resp1, resp2=resp2, labels=ptid_labels, title_str=title_str))
    fig_handles[0].set_size_inches(18, 11)
    fig_handles[1].set_size_inches(11, 11)
    return fig_handles

def plot_longitudinal_responses_by_ptid(resp_mat, prot_name=None, ptid=None, timepoint_labels=None, plot_type='line',
                                        subplot_flag=False, fig_size=(18, 11)):
    """
    plot responses of a single ptid from longitudinal timepoints
    
    Parameters:
    ----------
    resp_mat : np.ndarray, or pd.DataFrame
        antigen array responses where each row corresponds to responses for a single timepoint
    ptid : [ string | None (default)]
        ptid id (optional)
    timepoint_labels : [strings | None (default)]
        labels for each timepoint, used for legened plotting
    plot_type : ['line' (default) | 'bar']
        type of plot used to present data
    subplot_flag : [True | False]
        if True plot each timepoint in a different subplot, else plot all on a single plot.

    Returns:
    -------
    f : figure handle

    """
    N, num_antigens  = resp_mat.shape
   
    width = 0.8/N      # the width of the bars
    num_colors = max(3,N)
    cmap_name = "palettable.colorbrewer.qualitative." + 'Set1' + "_" + str(num_colors) + ".mpl_colors"
    colors = eval(cmap_name)

    if isinstance(resp_mat,pd.DataFrame):
        resp_mat = resp_mat.as_matrix()

    if subplot_flag:
        num_plots = N
    else:
        num_plots = 1
    f, axarr = plt.subplots(num_plots,1)
    f.set_tight_layout(True)
    f.set_size_inches(fig_size)

    for i in np.arange(N):
        
        if subplot_flag:
            curr_ax = axarr[i]
        else:
            curr_ax = axarr

        if plot_type == 'line':    
            curr_ax.plot(np.arange(1,num_antigens+1),resp_mat[i],color=colors[i])
        elif plot_type == 'bar':
            curr_ax.bar(np.arange(1,num_antigens+1), resp_mat[i], axis=1, width=width, color=colors[0])            
    
        if ptid is None:
            ptid = ""
        if prot_name is None:
            prot_name = ""
        curr_ax.set_title(prot_name + " longitudinal response for ptid " + ptid)
        curr_ax.set_yticks([])
        curr_ax.set_ylim(0, 60000)
        if subplot_flag:
            plt.legend(timepoint_labels[i])
    
    if not subplot_flag:
        plt.legend(timepoint_labels)

    return f


def plot_responses_by_clusters(data, cluster_labels, fig_size=(18,11), cluster_title_prefix=""):
    """ 
    Plot antigen array responses of individual samples by cluster assignment. 
    Each cluster is plotted in a sepearate panel (subplot)

    Parameters:
    ----------
    data : [np.ndarray | pd.DataFrame] 
        Data for each sample - each row is the responses of a single subject to a set of antigens on the array
    cluster_labels : np.ndarray of integers
        Array of integers with the cluster assignment of each group. 
        Clusters are assumed to be labeled from 1:num_clusters.
    fig_size : 2-tuple of integers default set to (18, 11)
        (width, height) of figure
    cluster_title_prefix : [str | None (default)]
        prefix to the cluster titles,  e.g. which protein antigens were used for clustering.

    Returns:
    -------
    f : figure handle

    """
    num_clusters = len(np.unique(cluster_labels))

    f, axarr = plt.subplots(num_clusters,1)    
    f.set_tight_layout(True)
    f.set_size_inches(fig_size)
    
    antigen_inds = np.arange(1,data.shape[1]+1)

    # plot clusters  
    for i in np.arange(1,num_clusters+1):
    
        curr_inds = [cluster_labels == i]
        curr_data = data[curr_inds][:]

        axarr[i].plot(np.arange(antigen_inds, curr_data)) 
        axarr[i].set_title(title_prefix + " cluster " + str(i+1) + " (n = " + str(len(curr_inds)) + ")")
        axarr[i].set_yticks([])

    return f

def plot_responses_by_groups(data, group_labels, groups=None, group_names=None, fig_size=(18,11), title_prefix=''):
        """ 
        Plot antigen array responses of individual samples by group assignment. 
        Each group is plotted in a different panel (subplot)

        Parameters:
        ----------
        data : [np.ndarray |pd.DataFrame] 
            Data - each row is the responses of a single subject to a set of antigens on the array.  
        group_labels : [np.array | pd.Series]
            vector with the group assignemt of each sample. Can be integer based or string based labels.        
        groups : [list | pd.Series | None (default)] 
            Names of groups to be plotted. If None, will plot all groups that exist in group_labels
        group_names : [strs | None (default)]
            Names of groups used for plot titles. If none values from 'groups' wil be used instead 
        fig_size : 2-tuple of integers default set to (18, 11)
            (width, height) of figure
        title_prefix : str
            Prefix of title used for group title

        Returns:
        -------
        f : figure handle

        """
        if groups is None:
            groups = np.unique(group_labels)
        if group_names is None:
            group_names = groups

        num_plots = len(groups)

        f, axarr = plt.subplots(num_plots,1)    
        f.set_tight_layout(True)
        f.set_size_inches(fig_size)
        
        antigen_inds = np.arange(1,data.shape[1]+1)
        # plot groups:
        for i, group in enumerate(groups):
        
            g_inds = mlab.find(group_labels == group)
            curr_data = data.iloc[g_inds]
            
            axarr[i].plot(antigen_inds, curr_data.T.as_matrix())
            axarr[i].set_title(title_prefix + group_names[i]  + " (n = " + str(len(g_inds)) + ")")
            axarr[i].set_yticks([])

        return f

def plot_clustering_dendrograms(Z_struct, prot_names, labels, fig_path=None, orientation='left', fig_size=(18,11), fig_prefix=None):
    """
    plot all clustering dendrograms using the specific set of proteins in prot_names
    
    Parameters:
    ----------
    Z_struct: dictionary
        clustering structure (matlab style) returned by linkage indexed by ind_dict.keys()
    prot_names: list
        list of strings of protein antigens from which peptides are on the array_data_filename
    labels: list
        list of labels of datapoints to label the dendrogram. 
    fig_path: [string | None (default)]
        if set, will save figures in fig_path. if None (default) figures are not saved.
    fig_size : 2-tuple of integers default set to (18, 11)
        (width, height) of figure
    """
    for p in prot_names:
        f, axarr = plt.subplots(1)
        f.set_size_inches(fig_size)        
        sch.dendrogram(Z=Z_struct[p], labels=labels, orientation=orientation)
        axarr.set_title(p)

        if fig_path is not None:
            f.set_tight_layout(True)
            if fig_prefix is None:
                fig_prefix = ''
            filename = "".join([fig_path, fig_prefix, p, "_dendrograms.png"])
            f.savefig(filename, dpi=200)

def plot_summary_stat_boxplots_by_clusters(arr_df, clusters, prot_names, arr_summary_stats, sample_inds=None, fig_prefix=None, fig_path=None, fig_size=(11,11)):
    """ 
    Plot boxplots of summary stats by clusters
    
    Parameters:
    ----------
    arr_df: pandas.DataFrame
        dataframe of array data
    clusters: dictionary
        cluster assignment of each datapoint indexed by ind_dict.keys()
    prot_names: list
        list of strings of protein antigens from which peptides are on the array_data_filename
    arr_summary_stats: list
        list of array summary statistics for which boxplots are to be generated.
    sample_inds: [bool array | None]
        boolean array specifying which samples to use. If None (default) will plot all.
    fig_path: [string | None (default)]
        if set, will save figures in fig_path. if None (default) figures are not saved.
    fig_size : 2-tuple of integers default set to (11, 11)
            (width, height) of figure   
    """
    if sample_inds is None:
        sample_inds = [True]*arr_df.shape[0]
    
    for p in prot_names:
        for assay in arr_summary_stats:
            f = plt.figure()    
            f.set_size_inches(fig_size)       
            mbp.myboxplot_by_labels(arr_df[sample_inds][assay], clusters[p])
            plt.title("".join([p, " clusters ", assay]))
            plt.xlabel('Cluster #')
            num_clusters = len(np.unique(clusters[p]))
            # save to file only if save_flag is on
            if fig_path is not None:
                filename = "".join([fig_path, fig_prefix, p, "_", assay, "_boxplots_by_clusters_n_", str(num_clusters), ".png"])
                f.savefig(filename, dpi=200)


def plot_summary_stat_boxplots_by_exp_groups(arr_df, arr_summary_stats, sample_inds=None, fig_path=None, fig_prefix=None, fig_size=(11,11)):
    """"
    Plot boxplots of summary stats by experimental groups
    
    Parameters:
    ----------
    arr_df: pandas.DataFrame
        dataframe of array data
    arr_summary_stats: list
        list of array summary statistics for which boxplots are to be generated.
    sample_inds: [bool array | None (default)]
        defines which samples to use in boxplots (allows selecting subset). 
        If set to None (default) will use all samples.
    fig_path: [string | None (default)]
        if set, will save figures in fig_path. if None (default) figures are not saved.   
    fig_prefix: [string | None]
        if a subset of the data is plotted, allows specifying a figure name prefix
        for saving to files. 
    fig_size : 2-tuple of integers default set to (11, 11)
            (width, height) of figure
    """ 
    if sample_inds is None:
        sample_inds = [True]*arr_df.shape[0]
    if fig_prefix is None:
        fig_prefix = ''

    for assay in arr_summary_stats:        
        f = plt.figure()
        f.set_size_inches(fig_size)
        mbp.myboxplot_by_labels(arr_df[sample_inds][assay], arr_df[sample_inds]['group'])
        plt.title("".join([fig_prefix, " ", assay, "_responses"]))
        plt.xlabel('group #')
    
        # save to file only if save_flag is on
        if fig_path is not None:
            f.set_tight_layout(True)
            filename = "".join([fig_path, fig_prefix, assay, "_boxplots_by_groups.png"])
            f.savefig(filename, dpi=200)

def plot_responses_by_exp_groups(arr_df, antigen_inds, exp_groups, fig_path=None, fig_prefix=None, fig_size=(18,11), y_lims=[0, 60000]):
    """
    plot responses by groups
    """ 
    num_groups = len(exp_groups)
    f, axarr = plt.subplots(num_groups,1)
    f.set_size_inches(fig_size)

    # plot groups
    if len(exp_groups) == 1:
        i=0
        axarr.plot(np.arange(len(antigen_inds)), arr_df[antigen_inds].iloc[np.where(arr_df['group'] == exp_groups[i])].T)
        axarr.set_title(fig_prefix + " " + exp_groups[i] + " (n = " + str(len(np.where(arr_df['group'] == exp_groups[i])[0])) + ")")
        axarr.set_yticks([])
        axarr.set_ylim(y_lims)
    else:
        for i in np.arange(num_groups):

            axarr[i].plot(np.arange(len(antigen_inds)), arr_df[antigen_inds].iloc[np.where(arr_df['group'] == exp_groups[i])].T)
            axarr[i].set_title(fig_prefix + " " + exp_groups[i] + " (n = " + str(len(np.where(arr_df['group'] == exp_groups[i])[0])) + ")")
            axarr[i].set_yticks([])
            axarr[i].set_ylim(y_lims)

    # save to file only if save_flag is on
    if fig_path is not None:
        f.set_tight_layout(True)
        filename = "".join([fig_path, fig_prefix,  "_raw_responses_by_groups.png"])
        f.savefig(filename, dpi=200)


def plot_median_responses_by_exp_groups(group_medians, exp_groups, prot_str, fig_path=None, fig_prefix=None, fig_size=(18,11), y_lims=[0, 30000]):
    """ plot bar plots of median responses to a given set of antigens by experimental group"""
    width = 0.3 # the width of the bars 
    colors = ['r','b', 'k', 'g', 'm']
            
    f, axarr = plt.subplots(len(exp_groups),1)
    f.set_size_inches(fig_size)
    
    if len(exp_groups) == 1:
        rects1 = axarr.bar(np.arange(len(group_medians[(exp_groups[0], prot_str)])), group_medians[(exp_groups[0], prot_str)], width=width, color=colors[0])
        axarr.set_title(exp_groups[0] + " " + prot_str)
        axarr.set_ylim(y_lims)
    else:
        for i in np.arange(len(exp_groups)):
            rects1 = axarr[i].bar(np.arange(len(group_medians[(exp_groups[i], prot_str)])), group_medians[(exp_groups[i], prot_str)], width=width, color=colors[i])
            axarr[i].set_title(exp_groups[i] + " " + prot_str)
            axarr[i].set_ylim(y_lims)


    if fig_path is not None:    
        f.set_tight_layout(True)
        if fig_prefix is None:
            fig_prefix =''
        filename = "".join([fig_path, fig_prefix, prot_str,  "_median_responses_by_treatment_group.png"])
        f.savefig(filename, dpi=200)


def plot_median_responses_by_clusters(arr_df, antigen_inds, num_clusters, clusters, fig_path=None, fig_prefix='', fig_size=(18,11), y_lims=[0, 60000]):
    """
    plot bar plot of median responses by clusters.
    """ 
    f, axarr = plt.subplots(num_clusters,1)
    f.set_size_inches(fig_size)

    # plot groups
    for i in np.arange(num_clusters):
        axarr[i].bar(np.arange(len(antigen_inds)), np.median(arr_df[antigen_inds].loc[clusters == i+1].T, axis=1))
        axarr[i].set_title(fig_prefix + " cluster " + str(i+1) + " (n = " + str(len(np.where([clusters == i+1])[0])) + ")")
        axarr[i].set_yticks(y_lims)
        axarr[i].set_ylim(y_lims)

    # save to file only if save_flag is on
    if fig_path is not None:
        f.set_tight_layout(True)
        filename = "".join([fig_path, fig_prefix,  "_median_responses_by_clusters_n_", str(num_clusters), ".png"])
        f.savefig(filename, dpi=200)

def plot_raw_responses_by_clusters(arr_df, antigen_inds, num_clusters, clusters, fig_path=None, fig_prefix='', fig_size=(18,11), y_lims=[0, 60000]):
    """
    plot bar plot of median responses by clusters.
    """ 
    f, axarr = plt.subplots(num_clusters,1)
    f.set_size_inches(fig_size)

    # plot groups
    for i in np.arange(num_clusters):
        axarr[i].plot(np.arange(len(antigen_inds)), arr_df[antigen_inds].loc[clusters == i+1].T)
        axarr[i].set_title(fig_prefix + " cluster " + str(i+1) + " (n = " + str(len(np.where([clusters == i+1])[0])) + ")")
        axarr[i].set_yticks([])
        axarr[i].set_ylim(y_lims)

    # save to file only if save_flag is on
    if fig_path is not None:
        f.set_tight_layout(True)
        filename = "".join([fig_path, fig_prefix,  "_median_responses_by_clusters_n_", str(num_clusters), ".png"])
        f.savefig(filename, dpi=200)
