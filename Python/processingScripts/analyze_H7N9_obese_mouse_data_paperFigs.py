
from __future__ import division, print_function
from scipy import io
import pandas as pd
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
from sklearn import metrics
import amplotlib as amp
import amutils as amutils
import fdr as fdr


ROOT_PATH = '/Users/thertz/Dropbox/HertzLab/'
#ROOT_PATH = 'C:/Users/tomer.hertz/Dropbox/HertzLab/'


SAVE_PATH = ROOT_PATH + 'ArrayData/Influenza/H7N9/'
FIG_PATH = ROOT_PATH + 'ArrayData/Influenza/H7N9/paperFigs/'

prot_names = ['SHA_ha', 'SHA_na', 'Cal_ha', 'Cal_na']
prot_strs = ['H7', 'N9', 'H1', 'N1']

adjuvants = ['PBS', 'Vac', 'AS03']  # 'MF59' dropped out here. types of adjuvants used on both Obest and WT mice
exp_group_prefixes = ['WT_pre_', 'Ob_pre_', 'WT_post_', 'Ob_post_']

group_labels = {}
l = 1
for e in exp_group_prefixes:
    for a in adjuvants:
        group_labels[e + a] = l
        l = l + 1


assays = ['HAI_H7', 'MN_H7', 'HAI_H1', 'MN_H1']
# arr_summary_stats = ['H7_mag', 'N9_mag', 'H1_mag', 'N1_mag']
arr_summary_stats = ['H7_mag', 'N9_mag', 'H7_breadth', 'N9_breadth']


timepoints = ['Pre', 'Post']

# read in the data on HAI MN etc. from Erik:
serologyDf = pd.read_table(ROOT_PATH + '/ArrayData/Influenza/H7N9/OBVAC_H7_PostBoost_Serology_Data_for_Tomer_v2.txt', header=1)
immunologyDf = pd.read_table(ROOT_PATH + '/ArrayData/Influenza/H7N9/OBVAC_H7_PostBoost_IgG_IgM_Data_for_Tomer_v2.txt', header=0)

# generate new index column to match the ptid values in the arrayData
serologyDf['ptid'] = serologyDf['Cage'].map(str) + serologyDf['Mouse'].map(str)
serologyDf.set_index('ptid', inplace=True)

immunologyDf['ptid'] = immunologyDf['Cage'].map(str) + immunologyDf['Mouse'].map(str)
immunologyDf.set_index('ptid', inplace=True)

# read in all mat files of array data, and strip them out from the matstructs into a dataframe:
exp_dates = {'08_22_2014', '08_25_2014', '09_04_2014', '09_05_2014'}
arr_df = None
for dni, dirName in enumerate(exp_dates):

    filename_list = glob.glob(ROOT_PATH + '/ArrayData/Influenza/H7N9/' + dirName + '/*.mat')
    print(dni, filename_list[0])
    d = io.loadmat(filename_list[0], struct_as_record=False)['arrayData']
    matstruct = d[0][0]
    ptids = [matstruct.ptids[0][i][0] for i in np.arange(matstruct.ptids[0].shape[0])]
    group_names = [matstruct.groupNames[0][i][0] for i in np.arange(matstruct.ptids[0].shape[0])]
    res = matstruct.responseMatrix[0][0]
    antigens = [matstruct.antigenNames[i][0][0] for i in np.arange(matstruct.antigenNames.shape[0])]

    columns = antigens + ['group', 'date']
    data = np.column_stack((res, np.asarray(group_names)))
    data = np.column_stack((data, np.asarray([dirName]*len(group_names))))

    if arr_df is None:
        arr_df = pd.DataFrame(data, index=ptids, columns=columns)
    else:
        arr_df = pd.concat((arr_df, pd.DataFrame(data, index=ptids, columns=columns)), axis=0)


arr_df['group_label'] = arr_df['group'].map(group_labels)
# add HAI, NT data into the dataframe
arr_df['HAI_H7'] = np.nan
arr_df['MN_H7'] = np.nan
arr_df['HAI_H1'] = np.nan
arr_df['MN_H1'] = np.nan

for t in ['pre', 'post']:
    arr_inds = arr_df.group.str.contains(t)
    serology_inds = arr_df.index[arr_inds].map(lambda s: s[1:])

    arr_df.loc[arr_inds, 'HAI_H7'] = serologyDf.loc[serology_inds]['HAI H7 ' + t.capitalize() + ' Boost'].values
    arr_df.loc[arr_inds, 'MN_H7'] = serologyDf.loc[serology_inds]['MN H7 ' + t.capitalize() + ' Boost'].values
    arr_df.loc[arr_inds, 'HAI_H1'] = serologyDf.loc[serology_inds]['HAI H1 ' + t.capitalize() + ' Boost'].values
    arr_df.loc[arr_inds, 'MN_H1'] = serologyDf.loc[serology_inds]['MN H1 ' + t.capitalize() + ' Boost'].values

# create dictionay of all exp groups in the data where data is indices into the arr_df
group_inds = {}
for e in exp_group_prefixes:
    for a in adjuvants:
        group_inds[e + a] = arr_df.loc[arr_df['group'] == e + a].index

# read FIM12 peptide data
pep_df = pd.read_table(ROOT_PATH + 'ArrayData/Influenza/FIM12/docs/FIM12_Strains_PeptidesPrinted_May2014.txt')
pep_df = pep_df.set_index('print_name')


# background subtraction: subtract maximal BSA response for each antigen
bg_responses = np.asarray(arr_df.loc['BSA'][antigens].max())
fg_inds = arr_df.index != 'BSA'
bg_df = arr_df.loc['BSA']
arr_df = arr_df[fg_inds]

curr_data = np.array(arr_df.as_matrix(columns=[antigens]), dtype=float) - bg_responses.T
curr_data[np.where(curr_data < 0)] = 0
arr_df[antigens] = curr_data
arr_df[antigens] = arr_df[antigens].astype(float)


# now subtract out the PBS response from the vaccine group responses to remove baseline repertoire responses
wt_baseline_responses = np.asarray(arr_df[arr_df.group == 'WT_post_PBS'][antigens].median())
ob_baseline_responses = np.asarray(arr_df[arr_df.group == 'Ob_post_PBS'][antigens].median())

# wt_baseline_responses = np.asarray(arr_df[arr_df.group == 'WT_pre_PBS'][antigens].median())
# ob_baseline_responses = np.asarray(arr_df[arr_df.group == 'Ob_pre_PBS'][antigens].median())

wt_rows = arr_df.group.str.startswith('WT')
ob_rows = arr_df.group.str.startswith('Ob')

curr_data_wt = np.array(arr_df[wt_rows].as_matrix(columns=[antigens]), dtype=float) - \
    wt_baseline_responses.T
curr_data_ob = np.array(arr_df[ob_rows].as_matrix(columns=[antigens]), dtype=float) - \
    ob_baseline_responses.T

curr_data_wt[np.where(curr_data_wt < 0)] = 0
curr_data_ob[np.where(curr_data_ob < 0)] = 0

arr_df.loc[wt_rows, antigens] = curr_data_wt
arr_df.loc[ob_rows, antigens] = curr_data_ob
arr_df[antigens] = arr_df[antigens].astype(float)

ind_dict = {}  # indices of cols for each strain:
for p, s in zip(prot_names, prot_strs):
    ind_dict[p] = [col for col in arr_df.columns if (col[:6] == p and col[7:].isdigit())]  # define indices into columns of all specific HA and NA seqs:
    arr_df[s + '_mag'] = arr_df[ind_dict[p]].sum(axis=1)  # insert new columns into dataframe for overall magnitude for each strain

    # breadth - response is positive if it is above the mean response of that antigen across all samples
    for col in ind_dict[p]:
        arr_df[col + '_binarized'] = arr_df[col].map(lambda s: 1 if s > 2000 else 0)

    arr_df[s + '_breadth'] = arr_df[[col for col in arr_df.columns if (col[:6] == p and col.endswith('_binarized'))]].sum(axis=1)


# compute some stats - correlations between overall magnitude and breadth and HAI and NT measurements:
time_dict = {}
time_dict['Pre'] = arr_df.group.str.contains('pre')
time_dict['Post'] = arr_df.group.str.contains('post')


# compute median responses by treatment group:
group_medians = {}
group_stds = {}
for p, s in zip(prot_names, prot_strs):
    for g in group_inds.keys():
        curr_df = arr_df.loc[group_inds[g]][ind_dict[p]]
        group_medians[g, s] = curr_df.apply(np.median, axis=0)
        group_stds[g, s] = curr_df.apply(np.std, axis=1)

# magnitude and breadth correlations with alternate assays.
for t in time_dict.keys():
    for assay in assays:

        r, p = scipy.stats.spearmanr(arr_df.loc[time_dict[t], 'H7_mag'], arr_df.loc[time_dict[t], assay])
        print("spearman corr between Shanghai HA magnitude and " + assay + " " + t + " vaccination r = {:.2f}, p = {:.4g}".format(r, p))

        r, p = scipy.stats.spearmanr(arr_df.loc[time_dict[t], 'H7_breadth'], arr_df.loc[time_dict[t], assay])
        print("spearman corr between Shanghai HA breadth and " + assay + " " + t + " vaccination r = {:.2f}, p = {:.4g}\n".format(r, p))

        # plot linear correlation plots
        f = amp.plot_linear_corr(arr_df.loc[time_dict[t], 'H7_mag'], arr_df.loc[time_dict[t], assay], x_label="H7_mag ",
                                 y_label=assay, title_str=t + " boost")
        filename = "".join([FIG_PATH, "H7_mag_vs_", assay, "_", t, "_boost_correlation.png"])
        f.savefig(filename, dpi=200)
        
        filename = "".join([FIG_PATH, "H7_mag_vs_", assay, "_", t, "_boost_correlation.eps"])
        f.savefig(filename, dpi=1000)



        f = amp.plot_linear_corr(arr_df.loc[time_dict[t], 'H7_breadth'], arr_df.loc[time_dict[t], assay], x_label="H7_breadth ",
                                 y_label=assay, title_str=t + " boost")
        filename = "".join([FIG_PATH, "H7_breadth_vs_", assay, "_", t, "_boost_correlation.png"])
        f.savefig(filename, dpi=200)

        filename = "".join([FIG_PATH, "H7_breadth_vs_", assay, "_", t, "_boost_correlation.eps"])   
        f.savefig(filename, dpi=1000)


#  cluster using Andrew's package:
num_clusters = 4
dMat = {}  # distance matrices
Z_struct = {}  # clustering struct
dend = {}  # dendrogram struct
clusters = {}  # cluster labels

# cluster obese and WT by adjuvant:
for a in adjuvants:
    curr_inds = group_inds['Ob_post_' + a].append(group_inds['WT_post_' + a])
    dMat[a] = {}
    Z_struct[a] = {}
    dend[a] = {}
    clusters[a] = {}
    for p in ['SHA_ha', 'SHA_na']:
        # use Andrew's package which allows clustering using Spearman distances (sch.linkage, and pdist do not support this for some reason, unlike Matlab)
        (dMat[a][p], Z_struct[a][p], dend[a][p]) = hcp.computeHCluster(arr_df.loc[curr_inds][ind_dict[p]], method='complete', metric='spearman')
        clusters[a][p] = sch.fcluster(Z_struct[a][p], t=num_clusters, criterion='maxclust')


# compute ranksum p-values for comparisons of HAI, microneut assays for the Shanghai and Cal strains comparing different treatment groups across the same adjuvant
p_values = {assay: np.zeros(shape=(len(adjuvants), len(exp_group_prefixes))) for assay in assays}
q_values = {assay: np.zeros(shape=(len(adjuvants), len(exp_group_prefixes))) for assay in assays}
stats_df = pd.DataFrame()

for assay in assays + arr_summary_stats:
    res = []
    ind_labels = []
    for t in ['pre', 'post']:
        for ad in adjuvants:
            res.append(scipy.stats.ranksums(arr_df.loc[group_inds['Ob_' + t + '_' + ad]][assay], arr_df.loc[group_inds['WT_' + t + '_' + ad]][assay]))
            ind_labels.append(t + '_' + ad)
    if(stats_df.empty):
        stats_df = pd.DataFrame(index=ind_labels)
    z_vals, p_vals = zip(*res)
    reject, q_vals, temp1, temp2 = sm.stats.multipletests(p_vals, method='fdr_bh', alpha=0.2)
    stats_df[assay + '_p_vals'] = p_vals
    stats_df[assay + '_q_vals'] = q_vals


# reorder columns by sorted order:
sorted_inds = np.argsort(stats_df.columns)
stats_df = stats_df[stats_df.columns[sorted_inds]]
stats_df.filter(regex='p_vals').to_csv(path_or_buf=SAVE_PATH + "Ob_vs_WT_stats_by_group.csv")


# analyze single antigens:

# filter antigens based on percent responders:
filter_dict = amutils.compute_label_blinded_antigen_filters(arr_df, ind_dict, filter_threshold=0.4, pos_threshold=2000)
sig_antigens = amutils.compare_single_antigens_by_groups(arr_df, ['Ob_post_AS03', 'WT_post_AS03'], filter_dict, pos_threshold=2000)

# now generate dataframes for sig antigens from proteins of interest:
writer = pd.ExcelWriter(FIG_PATH + 'sig_antigen_stats.xlsx')
for p, p_name in zip(prot_names, prot_strs):
    sig_antigens_df = pd.DataFrame(sig_antigens[p], columns=['Antigen', 'Fisher Counts', 'p-value', 'adjusted p-value', 'q-value'])
    sig_antigens_df.to_excel(writer, sheet_name=p_name)
writer.save()


#  Figure 1: PBS responses of the groups are different
#  Plot PBS responses of the two groups in bar plot
fontsize = 24
width = 0.3  # the width of the bars
colors = ['r', 'b']
for p, s in zip(prot_names, prot_strs):

    curr_inds = np.arange(1, len(ind_dict[p]) + 1)
    for t in ['Post']:  # time_dict.keys():
        f, axarr = plt.subplots(1, 1)
        f.set_tight_layout(True)

        # plot PBS responses by group

        # here we have to fish our indices from original matlab antigens list since we are
        # working with the matrix form of these and not the dataframe
        prot_inds = [i for i, element in enumerate(antigens) if element.startswith(p)]

        rects1 = axarr.bar(curr_inds, ob_baseline_responses[prot_inds], width=width, color=colors[0])
        rects2 = axarr.bar(curr_inds + width, wt_baseline_responses[prot_inds], width=width, color=colors[1])
        axarr.legend(["Ob", "WT"], fontsize=fontsize)
        axarr.set_title("PBS post boost " + s + " Responses", fontsize=fontsize)
        axarr.set_yticks([0, 10000, 20000])
        axarr.set_ylim([0, 25000])
        axarr.set_xlabel('peptide #', fontsize=fontsize)
        axarr.set_ylabel('MFI', fontsize=fontsize)
        axarr.tick_params(axis='both', which='major', labelsize=20)

    filename = "".join([FIG_PATH, s, "_PBS_responses_by_treatment_group.png"])
    f.set_size_inches(18, 11)
    f.savefig(filename, dpi=200)


    filename = "".join([FIG_PATH, s, "_PBS_responses_by_treatment_group.eps"])
    f.savefig(filename, dpi=1000)


# Plot figures responses by treatment group
# Figure 2
width = 0.3  # the width of the bars
colors = ['r', 'b']
fig_adjuvants = ['Vac', 'AS03']
for p, s in zip(prot_names, prot_strs):

    curr_inds = np.arange(1, len(ind_dict[p]) + 1)
    for t in ['Post']:  # time_dict.keys():
        f, axarr = plt.subplots(len(fig_adjuvants), 1)
        f.set_tight_layout(True)

        # plot groups
        for i in np.arange(len(fig_adjuvants)):

            Ob_name = "Ob_" + t.lower() + "_" + fig_adjuvants[i]
            rects1 = axarr[i].bar(curr_inds, group_medians[Ob_name, s], width=width, color=colors[0])

            WT_name = "WT_" + t.lower() + "_" + fig_adjuvants[i]
            rects2 = axarr[i].bar(curr_inds + width, group_medians[WT_name, s], width=width, color=colors[1])
            axarr[i].legend(["Ob", "WT"])
            axarr[i].set_title(fig_adjuvants[i] + " " + t + " Boost " + s + " Responses", fontsize=fontsize)
            axarr[i].set_yticks([0, 10000, 20000])
            axarr[i].set_ylim([0, 25000])
            if (i == len(fig_adjuvants) - 1):
                axarr[i].set_xlabel('peptide #', fontsize=fontsize)
            axarr[i].set_ylabel('MFI', fontsize=fontsize)
            axarr[i].tick_params(axis='both', which='major', labelsize=20)

    filename = "".join([FIG_PATH, s, "_responses_by_treatment_group.png"])
    f.set_size_inches(18, 11)
    f.savefig(filename, dpi=200)

    filename = "".join([FIG_PATH, s, "_responses_by_treatment_group.eps"])  
    f.savefig(filename, dpi=1000)


# Figure 3 - magnitude and breadth boxplots of all groups
# plot boxplots of all groups:


# only plot for PBS, Vac, and AS03:
curr_df = arr_df[(arr_df.group != 'WT_post_MF59') & (arr_df.group != 'Ob_post_MF59') &
                 (arr_df.group != 'WT_pre_MF59') & (arr_df.group != 'Ob_pre_MF59')]
curr_time_dict = {}
curr_time_dict['Pre'] = curr_df.group.str.contains('pre')
curr_time_dict['Post'] = curr_df.group.str.contains('post')


for assay in arr_summary_stats:
    for t in ['Post']:  # time_dict.keys():
        f, axarr = plt.subplots(1)
        f.set_tight_layout(True)
        mbp.myboxplot_by_labels(curr_df[curr_time_dict[t]][assay], curr_df[curr_time_dict[t]]['group'])
        axarr.set_title("".join([t, " ", assay.replace('_', ' '), " responses"]), fontsize=16)
        axarr.tick_params(axis='both', which='major', labelsize=14)
        axarr.set_yscale('log')

        filename = "".join([FIG_PATH, t, "_", assay, "_boxplots_by_groups.png"])
        f.savefig(filename, dpi=200)

        filename = "".join([FIG_PATH, t, "_", assay, "_boxplots_by_groups.eps"])    
        f.savefig(filename, dpi=1000)


# Figure 4, 5 and 6 - clustering dendrograms, median responses and summary stats of WT vs. Obese for each group:
for a in ['Vac', 'AS03']:
    curr_inds = group_inds['Ob_post_' + a].append(group_inds['WT_post_' + a])
    amp.plot_clustering_dendrograms(Z_struct=Z_struct[a], prot_names=['SHA_ha'], labels=arr_df.loc[curr_inds].group, fig_prefix=a + '_', fig_path=FIG_PATH)
    amp.plot_median_responses_by_clusters(arr_df=arr_df.loc[curr_inds], antigen_inds=ind_dict['SHA_ha'], num_clusters=4,
                                          clusters=clusters[a]['SHA_ha'], y_lims=[0, 25000], fig_prefix=a, fig_path=FIG_PATH)
    #amp.plot_raw_responses_by_clusters(arr_df=arr_df.loc[curr_inds], antigen_inds=ind_dict['SHA_ha'], num_clusters=4, clusters=clusters[a]['SHA_ha'])
    amp.plot_summary_stat_boxplots_by_clusters(arr_df=arr_df.loc[curr_inds], clusters=clusters[a], prot_names=['SHA_ha'], 
                                               arr_summary_stats=['H7_breadth', 'H7_mag'], fig_prefix = a + '_', fig_path=FIG_PATH)
    


# # plot correlation matrix with dendrogram overlayed:
# for p in ['SHA_ha', 'SHA_na']:
#     f = plt.figure()
#     colInd = hcp.plotHColCluster(arr_df[post_inds][ind_dict[p]].T, method='complete', metric='spearman', titleStr=p, vRange=(0, 1),
#                                  col_labels=arr_df[post_inds].group)
#     f = plt.gcf()
#     filename = "".join([FIG_PATH + p + "postBoost_clustering_matrix_and_dendrogram"])
#     f.savefig(filename, dpi=200)
