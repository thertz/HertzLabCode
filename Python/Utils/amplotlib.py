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


def plot_responses_by_ptids(resp_mat, ptid_labels=None, plot_type='line', subplot_flag=False):
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
            curr_ax.plot(curr_inds, resp_mat[i],color=colors[i])
            
        elif plot_type == 'bar':
            curr_ax.bar(curr_inds, resp_mat[i], width=width, color=colors[i])            
    
        if (ptid_labels is not None) and subplot_flag:
            curr_ax.set_title("Response for ptid " + ptid_labels[i])

        #curr_ax.set_yticks([])
        curr_ax.set_ylim(0, 60000)
    
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

def plot_longitudinal_responses_by_ptid(resp_mat, ptid=None, timepoint_labels=None, plot_type='line',
                                        subplot_flag=False):
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
    f.set_size_inches(18,11)

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
        curr_ax.set_title("longitudinal response for ptid " + ptid)
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

def plot_responses_by_groups(data, group_labels, groups=None, group_names=None, fig_size=(18,11), title_prefix=""):
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
            
            axarr[i].plot(antigen_inds,curr_data.T.as_matrix())
            axarr[i].set_title(title_prefix + group_names[i]  + " (n = " + str(len(g_inds)) + ")")
            axarr[i].set_yticks([])

        return f




