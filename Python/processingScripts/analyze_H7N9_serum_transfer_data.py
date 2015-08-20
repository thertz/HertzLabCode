
from __future__ import division, print_function
from scipy import io
import numpy as np
import myboxplot as mbp
import scipy
import scipy.stats
import glob
import hclusterplot as hcp
import myboxplot as mbp
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import statsmodels.api as sm
from sklearn import metrics, cross_validation, linear_model, ensemble
import amplotlib as amp
import pandas as pd

ROOT_PATH = '/Users/thertz/Dropbox/HertzLab/'
SAVE_PATH = ROOT_PATH + 'ArrayData/Influenza/H7N9/'
FIG_PATH = ROOT_PATH + 'ArrayData/Influenza/H7N9/serumTransferFigsRepeat2/'
#plt.rcParams.update({'font.size': 12})

prot_names = ['SHA_ha', 'SHA_na', 'Cal_ha', 'Cal_na']
prot_strs = ['H7', 'N9', 'H1', 'N1']

adjuvants = ['']  # types of adjuvants used on both Obest and WT mice
exp_group_prefixes = ['Obese', 'Normal', 'Obese_pool', 'Normal_pool']

group_labels = {}
l = 1
for e in exp_group_prefixes:
    for a in adjuvants:
        group_labels[e + a] = l
        l = l + 1


# assays = ['HAI_H7', 'MN_H7',  'HAI_H1', 'MN_H1', ]
arr_summary_stats = ['H7_mag', 'N9_mag']#, 'H1_mag', 'N1_mag']

# read in mat file of array data, and strip them out from the matstructs into a dataframe:
arr_df = None
#fileName = '/Users/thertz/Dropbox/HertzLab/ArrayData/Influenza/H7N9/06_03_2015/H7N9_06_03_15_arrayData_median.mat'
#fileName = '/Users/thertz/Dropbox/HertzLab/ArrayData/Influenza/H7N9/06_10_2015/H7N9_06_10_15_arrayData_median_532.mat'
fileName = '/Users/thertz/Dropbox/HertzLab/ArrayData/Influenza/H7N9/06_18_2015/H7N9_06_18_15_arrayData_median.mat'

d = io.loadmat(fileName, struct_as_record=False)['arrayData']
matstruct = d[0][0]
ptids = [matstruct.ptids[0][i][0] for i in np.arange(matstruct.ptids[0].shape[0])]
group_names = [matstruct.groupNames[0][i][0] for i in np.arange(matstruct.ptids[0].shape[0])]
res = matstruct.responseMatrix[0][0]
antigens = [matstruct.antigenNames[i][0][0] for i in np.arange(matstruct.antigenNames.shape[0])]

columns = antigens + ['group']
data = np.column_stack((res, np.asarray(group_names)))
arr_df = pd.DataFrame(data, index=ptids, columns=columns)


# hack hack hack
#arr_df.loc['Obese_pool','group'] = 'Obese_pool'
#arr_df.loc['Normal_pool','group'] = 'Normal_pool'

# hack hack hack for repeat data 
#arr_df.loc['Obese_Pool_new','group'] = 'Obese_pool'
#arr_df.loc['Normal_Pool_new','group'] = 'Normal_pool'

arr_df.loc['Obese_Pool_AS','group'] = 'Obese_pool'
arr_df.loc['Normal_Pool_AS','group'] = 'Normal_pool'



arr_df['group_label'] = arr_df['group'].map(group_labels)


# create dictionay of all exp groups in the data where data is indices into the arr_df
group_inds = {}
for e in exp_group_prefixes:
    for a in adjuvants:
        group_inds[e + a] = arr_df.loc[arr_df['group'] == e + a].index

# read FIM12 peptide data
pep_df = pd.read_table(ROOT_PATH + 'ArrayData/Influenza/FIM12/docs/FIM12_Strains_PeptidesPrinted_May2014.txt')
pep_df = pep_df.set_index('print_name')


# background subtraction: subtract maximal BSA response for each antigen
# bg_responses = np.asarray(arr_df.loc['BSA'][antigens].max())
# fg_inds = arr_df.index != 'BSA'
# bg_df = arr_df.loc['BSA']
# arr_df = arr_df[fg_inds]

# bg_responses = np.asarray(arr_df.loc['BSA_new'][antigens].max())
# fg_inds = arr_df.index != 'BSA_new'
# bg_df = arr_df.loc['BSA_new']
# arr_df = arr_df[fg_inds]

bg_responses = np.asarray(arr_df.loc['BSA_AS'][antigens].max())
fg_inds = arr_df.index != 'BSA_AS'
bg_df = arr_df.loc['BSA_AS']
arr_df = arr_df[fg_inds]


curr_data = np.array(arr_df.as_matrix(columns=[antigens]), dtype=float) - bg_responses.T
curr_data[np.where(curr_data < 0)] = 0
arr_df.antigens = curr_data
arr_df[antigens] = arr_df[antigens].astype(float)

ind_dict = {}  # indices of cols for each strain:
for p, s in zip(prot_names, prot_strs):
    ind_dict[p] = [col for col in arr_df.columns if (col[:6] == p and col[7:].isdigit())]  # define indices into columns of all specific HA and NA seqs:
    arr_df[s + '_mag'] = arr_df[ind_dict[p]].sum(axis=1)  # insert new columns into dataframe for overall magnitude for each strain

    # breadth - response is positive if it is above the mean response of that antigen across all samples
    for col in ind_dict[p]:
        arr_df[col + '_binarized'] = arr_df[col].map(lambda s: 1 if s > 2000 else 0)

    arr_df[s + '_breadth'] = arr_df[[col for col in arr_df.columns if (col[:6] == p and col.endswith('_binarized'))]].sum(axis=1)

# compute median responses by treatment group:
group_medians = {}
group_stds = {}
for p, s in zip(prot_names, prot_strs):
    for g in group_inds.keys():
        curr_df = arr_df.loc[group_inds[g]][ind_dict[p]]
        group_medians[g, s] = curr_df.apply(np.median, axis=0)
        group_stds[g, s] = curr_df.apply(np.std, axis=1)

#  cluster using Andrew's package:
num_clusters = 4
dMat = {}  # distance matrices
Z_struct = {}  # clustering struct
dend = {}  # dendrogram struct
clusters = {}  # cluster labels
cluster_treatment_stats = {}
pred_treatment_labels = {}  # predicted treatment labels based on clustering

#Plot response per group
for p in ['SHA_ha', 'SHA_na']:
    f = amp.plot_responses_by_groups(arr_df[ind_dict[p]], arr_df.group)
    filename = "".join([FIG_PATH, p,  "_responses_by_groups.png"])
    f.savefig(filename, dpi=200)

# plot individual responses of obese transfer and normal transfer:
for p in ['SHA_ha', 'SHA_na']:
    for g in ['Obese', 'Normal']:
        f, axarr = plt.subplots(group_inds[g].shape[0], 1)
        f.set_tight_layout(True)
        for i, curr_g in enumerate(group_inds[g]):
            axarr[i].plot(arr_df.loc[curr_g][ind_dict[p]])
            axarr[i].set_ylim(0, 60000)
            axarr[i].set_yticks([])
            axarr[i].set_title(p + " " + curr_g)
        filename = "".join([FIG_PATH, p, "_", g, "_responses_by_ptids.png"])
        f.savefig(filename, dpi=200)

# cluster data:
for k in ind_dict.keys():
    # use Andrew's package which allows clustering using Spearman distances (sch.linkage, and pdist do not support this for some reason, unlike Matlab)
    (dMat[k], Z_struct[k], dend[k]) = hcp.computeHCluster(arr_df[ind_dict[k]], method='complete', metric='spearman')
    clusters[k] = sch.fcluster(Z_struct[k], t=num_clusters, criterion='maxclust')


# Plot figures for a given clustering solution - currently only performed for the Shanghai strain:
for p in ['SHA_ha', 'SHA_na']:
    f, axarr = plt.subplots(num_clusters, 1)
    f.set_tight_layout(True)
    f.set_size_inches(18, 11)
    # plot clusters
    for i in np.arange(num_clusters):

        axarr[i].plot(np.arange(len(ind_dict[p])), arr_df[ind_dict[p]].loc[clusters[p] == i+1].T)
        axarr[i].set_title(p + " cluster " + str(i+1) + " (n = " + str(len(np.where([clusters[p] == i+1])[0])) + ")")
        axarr[i].set_yticks([])
        axarr[i].legend(arr_df.loc[clusters[p] == i+1].index)
    filename = "".join([FIG_PATH, p,  "_responses_by_clusters_n_", str(num_clusters), ".png"])
    f.savefig(filename, dpi=200)

# plot cluster dendrogram:
for p in ['SHA_ha', 'SHA_na']:
    f, axarr = plt.subplots(1)
    f.set_tight_layout(True)
    f.set_size_inches(18,11)
    sch.dendrogram(Z_struct[p],color_threshold=np.inf, labels=arr_df.index, orientation='left')
    axarr.set_title(p)
    filename = "".join([FIG_PATH, p, "_dendrograms.png"])
    f.savefig(filename, dpi=200)


# plot boxplots of all clusters
for p in ['SHA_ha', 'SHA_na']:
    for sum_stat in arr_summary_stats:
        f = plt.figure()
        f.set_size_inches(18, 11)
        mbp.myboxplot_by_labels(arr_df[sum_stat], clusters[p])
        plt.title("".join([p, " clusters ", sum_stat]))
        plt.xlabel('Cluster #')
        filename = "".join([FIG_PATH, p, "_", sum_stat, "_boxplots_by_clusters_n_", str(num_clusters), ".png"])
        f.savefig(filename, dpi=200)

# # plot boxplots of all groups:
# for assay in assays + arr_summary_stats:
#     for t in time_dict.keys():
#         f = figure()
#         f.set_size_inches(18, 11)
#         mbp.myboxplot_by_labels(arr_df[time_dict[t]][assay], arr_df[time_dict[t]]['group'])
#         plt.title("".join([t, "_", assay, "_responses"]))
#         plt.xlabel('group #')
#         filename = "".join([FIG_PATH, t, "_", assay, "_boxplots_by_groups.png"])
#         f.savefig(filename, dpi=200)


# Plot figures responses by treatment group
# width = 0.3 # the width of the bars
# colors = ['r','b']
# for p, s in zip(prot_names, prot_strs):

#     curr_inds = np.arange(1,len(ind_dict[p])+1)
#     for t in time_dict.keys():
#         f, axarr = plt.subplots(len(adjuvants),1)
#         f.set_tight_layout(True)
#         f.set_size_inches(18,11)
#         # plot clusters
#         for i in np.arange(len(adjuvants)):

#             Ob_name = "Ob_" + t.lower() + "_" + adjuvants[i]
#             rects1 = axarr[i].bar(curr_inds, group_medians[Ob_name, s], width=width, color=colors[0])

#             WT_name = "WT_" + t.lower() + "_" + adjuvants[i]
#             rects2 = axarr[i].bar(curr_inds+width, group_medians[WT_name, s], width=width, color=colors[1])
#             axarr[i].legend(["Ob", "WT"])
#             axarr[i].set_title(adjuvants[i] + " " + t + " Boost " + s + " Responses")
#             axarr[i].set_yticks([])

#     filename = "".join([FIG_PATH, s,  "_responses_by_treatment_group.png"])
#     f.savefig(filename, dpi=200)

# # plot samples longitudinaly for each ptid:
#     ptids = arr_df.index.map(lambda s: s[1:])
# uniq_ptids = np.unique(ptids)
# for p in uniq_ptids:
#     ptid_inds = mlab.find(ptids == p)
#     if len(ptid_inds) == 1:
#         continue
#     resp_df = arr_df.iloc[ptid_inds][ind_dict['SHA_ha']]
#     timepoints = arr_df.iloc[ptid_inds]['group']
#     amp.plot_longitudinal_responses_by_ptid(resp_df, ptid=p, timepoint_labels=timepoints, plot_type='line', subplot_flag=False)


#Plot median response per cluster
# for p in ['SHA_ha', 'SHA_na']:
#     f, axarr = plt.subplots(num_clusters,1)
#     f.set_tight_layout(True)
#     f.set_size_inches(18,11)
#     # plot clusters
#     for i in np.arange(num_clusters):

#         axarr[i].bar(np.arange(len(ind_dict[p])), np.median(arr_df[post_inds][ind_dict[p]].loc[clusters[p] == i+1].T, axis=1))
#         axarr[i].set_title(p + " cluster " + str(i+1) + " (n = " + str(len(np.where([clusters[p] == i+1])[0])) + ")")
#         axarr[i].set_yticks([])
#         axarr[i].set_ylim(0, 20000)

#     filename = "".join([FIG_PATH, p,  "_median_responses_by_clusters_n_", str(num_clusters), ".png"])
#     f.savefig(filename, dpi=200)

# # Plot median response per pairs of clusters of interest
# for p in ['SHA_ha']:

#     ind = np.arange(len(ind_dict[p]))+1  # the x locations for the groups
#     width = 0.3       # the width of the bars

#     colors = ['r','b']
#  # plot pairs of clusters with significant differences
#     for i in sig_arr_df_HA.index:

#         f, ax = plt.subplots()
#         f.set_tight_layout(True)
#         f.set_size_inches(18,11)

#         #med_resp = np.median(arr_df[ind_dict[p]].loc[clusters[p] == i[0]].T, axis=1)
#         rects1 = ax.bar(ind, med_resp[:,i[0]-1], width, color=colors[0])
#         #med_resp = np.median(arr_df[ind_dict[p]].loc[clusters[p] == i[1]].T, axis=1)
#         rects2 = ax.bar(ind+width, med_resp[:,i[1]-1], width, color=colors[1])
#         ax.legend(i)
#         ax.set_title(p + " median responses of clusters")
   #         filename = "".join([FIG_PATH, p,  "_bar_median_responses_clusters_", str(i[0]), "_",str(i[1]), ".png"])
#         f.savefig(filename, dpi=200)


# plot correlation matrix of HA and NA clustering on Shanghai
# for p in ['SHA_ha', 'SHA_na']:

#     col_ind = np.argsort(clusters[p])
#     f = plt.figure()
#     ax_matrix = f.add_axes([0.3,0.1,0.6,0.6])
#     my_norm = mpl.colors.Normalize(vmin=0.4, vmax=1)
#     im = ax_matrix.matshow(dMat[p][colInd,:][:,colInd].T,aspect='auto', origin='upper', cmap=cm.RdBu_r,norm=my_norm)
#     ax_matrix.set_xticks([])

#     c_inds = []
#     c_inds = [len(np.where(clusters[p]==i)[0]) for i in np.arange(1,num_clusters+1)]
#     ax_matrix.set_yticks(np.cumsum(c_inds)[:-1])
#     ax_matrix.set_xticks(np.cumsum(c_inds)[:-1])
#     # Plot colorbar.
#     ax_color = f.add_axes([0.91,0.1,0.02,0.6])
#     plt.colorbar(im, cax=ax_color)


# plot correlation matrix with dendrogram overlayed:
for p in ['SHA_ha', 'SHA_na']:
    colInd = hcp.plotHColCluster(arr_df[ind_dict[p]].T,method='complete', metric='spearman',titleStr=p,vRange=(0, 1))
    #, \col_labels=arr_df[post_inds].group)
    f = plt.gcf()
    filename = "".join([FIG_PATH + p + "postBoost_clustering_matrix_and_dendrogram"])
    f.savefig(filename, dpi=200)
