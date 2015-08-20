
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


SAVE_PATH = ROOT_PATH + 'ArrayData/Influenza/FIM12/'
FIG_PATH  = ROOT_PATH + 'ArrayData/Influenza/FIM12/comparisonFigs/'
#plt.rcParams.update({'font.size': 12})


# read in the data on HAI MN etc. from Sanofi on FIM12:
rawDf = pd.read_csv(ROOT_PATH + '/ArrayData/Influenza/FIM12/FHCRC_6Nov14.csv')
renameD = {'U_HI_YNU': 'hi_dose',
           'U_SI_YNU': 'seasonal'}

rawDf = rawDf.rename_axis(renameD,axis=1)
rawDf['hi_dose'] = rawDf.hi_dose.map({'No':0.,'Yes':1.,'Unknown':2})
#rawDf['SUB_ID'] = rawDf.SUB_ID.map(lambda s: s[4:])
#rawDf['SUB_ID'] = rawDf.SUB_ID.str.slice(4)
rawDf = rawDf.set_index('SUB_ID') # Index is now not a column, but rather an index into the dataframe.

# print some stats:
print(rawDf.shape)
print (rawDf.columns)
print (rawDf.index[:10])


# read in all mat files of array data, and strip them out from the matstructs into a dataframe:
arrDf = None
for fni,fileName in enumerate(glob.glob(ROOT_PATH + '/ArrayData/Influenza/FIM12/*.mat')):
    print(fni, fileName)
    d = io.loadmat(fileName,struct_as_record=False)['arrayData']
    matstruct = d[0][0]
    ptids = [matstruct.ptids[0][i][0] + '_%d' % fni for i in np.arange(matstruct.ptids[0].shape[0])]
    res = matstruct.responseMatrix[0][0]
    antigens = [matstruct.antigenNames[i][0][0] for i in np.arange(matstruct.antigenNames.shape[0])]
    if arrDf is None:
        arrDf = pd.DataFrame(res,index=ptids,columns=antigens)
    else:
        arrDf = pd.concat((arrDf,pd.DataFrame(res,index=ptids,columns=antigens)),axis=0)

# reverse order of the ptid to allow joining with the other dataframe
arrDf['ptid'] = arrDf.index.map(lambda s: '-'.join(s.split('_')[:2][::-1]))
arrDf.set_index('ptid', inplace=True)
    
# read in case-control meta-data file sent with FIM12 serum samples
md_df  = pd.read_table(ROOT_PATH + 'ArrayData/Influenza/FIM12/docs/FIM12_mircoarray_samples_12Dec13.txt')
md_df.set_index('sub_id',inplace=True) # index is now not a column, but rather index into dataframe

# add label columns for case-control and for treatment assignment
md_df['case_control_label'] = md_df['case_control'].map(lambda s: 1 if s == 'Case' else 0)
md_df['treatment_label'] = md_df['trta'].map(lambda s: 1 if s == 'Fluzone High-Dose' else 0)

#arrDf.PTID.isin(rawDf.index).sum()
#arrDf.PTID.map(lambda s: s in rawDf.index)

#print arrDf.shape
#print arrDf.columns
#print arrDf.index[:10]

# inner join on two dataframes, allows us to now have one in which the index keys have all been merged!
df = arrDf.join(rawDf,how='inner')
df = df.join(md_df, how = 'inner')

# read FIM12 peptide data
pep_df = pd.read_table(ROOT_PATH + 'ArrayData/Influenza/FIM12/docs/FIM12_Strains_PeptidesPrinted_May2014.txt')
pep_df = pep_df.set_index('print_name')

prot_names = ['Vic_HA', 'Vic_NA','Cal_HA', 'Cal_NA', 'Wis_HA', 'Wis_NA']
ind_dict = {}; # indices of cols for each strain:
for p in prot_names:
    ind_dict[p] = [col for col in df.columns if (col[:6] == p and col[7:].isdigit())] # define indices into columns of all specific HA and NA seqs:
    df[p + '_mag']     = df[ind_dict[p]].sum(axis=1) # insert new columns into dataframe for overall magnitude for each strain

    # breadth - response is positive if it is above the mean response of that antigen across all samples
    for col in ind_dict[p]:
        df[col + '_binarized'] = df[col].map(lambda s: 1 if s > df[col].mean() else 0)

    df[p + '_breadth'] = df[[col for col in df.columns if (col[:6]== p and col.endswith('_binarized'))]].sum(axis=1)


# compute some stats - correlations between overall magnitude and breadth and HAI, NT and ELLA measurements:

# magnitude
r,p = scipy.stats.spearmanr(df.Vic_HA_mag,df.HAI_A_VICTORIA_361_2011_H3N2_CELL)
print("spearman corr between HA magnitude and HAI r = {:.2f}, p = {:.2f}".format(r,p))

r,p = scipy.stats.spearmanr(df.Vic_HA_mag,df.NT_A_VICTORIA_361_2011_H3N2_CELL)
print("spearman corr between HA magnitude and NT r = {:.2f}, p = {:.2f}".format(r,p))

r,p = scipy.stats.spearmanr(df.Vic_NA_mag,df.ELLA_N2)
print("spearman corr between NA magnitude and ELLA r = {:.2f}, p = {:.2f}".format(r,p))


# breadth:
r,p = scipy.stats.spearmanr(df.Vic_HA_breadth,df.HAI_A_VICTORIA_361_2011_H3N2_CELL)
print("spearman corr between HA magnitude and HAI r = {:.2f}, p = {:.2f}".format(r,p))

r,p = scipy.stats.spearmanr(df.Vic_HA_breadth,df.NT_A_VICTORIA_361_2011_H3N2_CELL)
print("spearman corr between HA magnitude and NT r = {:.2f}, p = {:.2f}".format(r,p))

r,p = scipy.stats.spearmanr(df.Vic_NA_breadth,df.ELLA_N2)
print("spearman corr between NA magnitude and ELLA r = {:.2f}, p = {:.2f}".format(r,p))


# now cluster the data using linkage: 

#Z_struct_Vic = cluster.hierarchy.linkage(df[[col for col in df.columns if col[:6] == 'Vic_HA']],method='complete',metric='correlation')
#Z_struct_Vic = cluster.hierarchy.linkage(df[[col for col in df.columns if col[:6] == 'Vic_NA']],method='complete',metric='correlation')

# cluster using Andrew's package:
num_clusters            = 6
dMat                    = {} # distance matrices
Z_struct                = {} # clustering struct
dend                    = {} # dendrogram struct
clusters                = {} # cluster labels
cluster_treatment_stats = {}
cluster_infection_stats = {}
pred_treatment_labels   = {} # predicted treatment labels based on clustering
pred_infection_labels   = {} # predicted infection labels based on clustering

for k in ind_dict.keys():
    # use Andrew's package which allows clustering using Spearman distances (sch.linkage, and pdist do not support this for some reason, unlike Matlab)
    (dMat[k], Z_struct[k], dend[k]) = hcp.computeHCluster(df[ind_dict[k]],method='complete',metric='spearman')
    clusters[k] = sch.fcluster(Z_struct[k], t = num_clusters, criterion='maxclust')
    
    # compute cluster homgeneity and completness (purity and accuracy) for treatment label and for infection status:
    pred_treatment_labels[k] = np.zeros(shape=(df.shape[0]))
    pred_infection_labels[k] = np.zeros(shape=(df.shape[0]))
    for i in np.arange(num_clusters):
        
        c_inds = np.where(clusters[k] == i)
        
        cluster_treatment_label = (1 if df['treatment_label'].values[c_inds].sum()/len(df['treatment_label'].values[c_inds]) >= 0.5 else 0)
        pred_treatment_labels[k][c_inds] = cluster_treatment_label

        cluster_infection_label = (1 if df['case_control_label'].values[c_inds].sum()/len(df['case_control_label'].values[c_inds]) >= 0.5 else 0)
        pred_infection_labels[k][c_inds] = cluster_infection_label

    cluster_treatment_stats[k] = metrics.homogeneity_completeness_v_measure(df['treatment_label'].values, pred_treatment_labels[k])
    cluster_infection_stats[k] = metrics.homogeneity_completeness_v_measure(df['case_control_label'].values, pred_treatment_labels[k])


# compute one way anova over clustering solutions:
for p in ['Vic_HA', 'Vic_NA']:    
    for assay in Vic_assays:
        group_samples = {}
        for i in arange(1,num_clusters+1):
            group_samples[i] = np.asarray(df[assay].loc[clusters[p] == i])
            group_samples[i] = group_samples[i][~np.isnan(group_samples[i])]
        (F, p_anova) = scipy.stats.f_oneway(*group_samples.values())
        print(p, assay, F, p_anova)
# compute ranksum p-values for comparisons of HAI, microneut and ELLA assays for the H3N2 strain comparing different clusters:
Vic_assays = ['HAI_A_VICTORIA_361_2011_H3N2_CELL', 'NT_A_VICTORIA_361_2011_H3N2_CELL', 'ELLA_N2' ]
assay_strs = {Vic_assays[0]:'HAI', Vic_assays[1]:'NT', Vic_assays[2]:'ELLA'}
prot_stats = {}
for p in ['Vic_HA', 'Vic_NA']:

    p_values = {assay: np.zeros(shape=(num_clusters, num_clusters)) for assay in Vic_assays}
    q_values = {assay: np.zeros(shape=(num_clusters, num_clusters)) for assay in Vic_assays}
    stats_df = pd.DataFrame()
    
    for assay in Vic_assays:
        res  = []
        c_inds = []
        for i in np.arange(num_clusters):
            for j in np.arange(i+1, num_clusters):
                res.append(scipy.stats.ranksums(df[assay].loc[clusters[p] == i+1], df[assay].loc[clusters[p] == j+1]))
                c_inds.append((i+1, j+1))
        stats_df['cluster_inds'] = c_inds
        stats_df = stats_df.set_index('cluster_inds')
        z_vals, p_vals = zip(*res)   
        reject, q_vals, temp1, temp2 = sm.stats.multipletests(p_vals,method='fdr_bh',alpha=0.2)        
        stats_df[assay_strs[assay] + '_p_vals'] = p_vals;
        stats_df[assay_strs[assay] + '_q_vals'] = q_vals;
        prot_stats[p] = stats_df
        
        # list comprehension syntax for the loop on clusters above, before we had multiple commands in the loop
        #res = [scipy.stats.ranksums(df[assay].loc[clusters[p] == i+1], df[assay].loc[clusters[p] == j+1]) \
        #    for i in np.arange(num_clusters) for j in np.arange(i+1, num_clusters)]
        
# identify all pairs of clusters with q-values < 0.2:
sig_df_HA = prot_stats['Vic_HA'].loc[(prot_stats['Vic_HA']['HAI_q_vals']<0.2) | (prot_stats['Vic_HA']['NT_q_vals']<0.2)]
sig_df_NA = prot_stats['Vic_NA'].loc[(prot_stats['Vic_NA']['ELLA_q_vals']<0.2)]

# collect median responses of each cluster to identify specific antigens in which responses in one group are greater than the other:
p = 'Vic_HA'
med_resp = []
for i in arange(num_clusters):
    med_resp.append(np.median(df[ind_dict[p]].loc[clusters[p] == i+1].T, axis=1))
med_resp = np.column_stack(np.asarray(med_resp))

# we are using the knowldege that the responses in cluster 6 are greater than 1 and 5 for HAI and NT accordingly here
for i in sig_df_HA.index: 
    if i == (5,6): # only 3-5 antigens here at all!
        d_inds = (np.where(med_resp[:,i[0]-1] < med_resp[:,i[1]-1]) and np.where(med_resp[:,i[1]-1] - med_resp[:,i[0]-1] > 1000))[0]
    elif i == (1,2):
        d_inds = (np.where(med_resp[:,i[0]-1] < med_resp[:,i[1]-1]) and np.where(med_resp[:,i[1]-1] - med_resp[:,i[0]-1] > 5000))[0]
    else:
        d_inds = (np.where(med_resp[:,i[0]-1] < med_resp[:,i[1]-1]) and np.where(med_resp[:,i[1]-1] - med_resp[:,i[0]-1] > 10000))[0]
    d_strs = []
    filename = SAVE_PATH + p + "_sig_peptide_list_clusters_" + str(i[0]) + "_" + str(i[1]) + ".txt"
    for d in d_inds:
        d_strs.append('Vic_HA_' + str(d+1))
    diff_df = pep_df[pep_df.index.isin(d_strs)].copy()
    diff_df.drop('modified', axis=1, inplace=True)
    print(diff_df)
    diff_df.to_csv(path_or_buf=filename,sep='\t')

# now translate indices into cluster inds: - not used since code modified to include tuple of indices into dataframe as index!
#c_inds = np.unravel_index(sig_df_HA.index,(num_clusters, num_clusters-1))


print("Vic HA clusters with significant differences in HAI or NT assays")
print(sig_df_HA)

# plot boxplots of all clusters
for p in ['Vic_HA', 'Vic_NA']:
    for assay in Vic_assays:
        f = figure()
        mbp.myboxplot_by_labels(df[assay], clusters[p])
        plt.title("".join([p, " clusters ", assay_strs[assay]]))
        plt.xlabel('Cluster #')
        filename = "".join([FIG_PATH, p, "_", assay_strs[assay], "_boxplots_by_clusters_n_", str(num_clusters), ".png"])
        f.savefig(filename, dpi=200)    
     
# Plot figures for a given clustering solution - currently only performed for the H3N2 victoria strain:
for p in ['Vic_HA', 'Vic_NA']:
    
    f, axarr = plt.subplots(num_clusters,1)    
    f.set_tight_layout(True)
    f.set_size_inches(18,11)
    # plot clusters  
    for i in np.arange(num_clusters):
        
        axarr[i].plot(np.arange(len(ind_dict[p])), df[ind_dict[p]].loc[clusters[p] == i+1].T)
        axarr[i].set_title(p + " cluster " + str(i+1) + " (n = " + str(len(np.where([clusters[p] == i+1])[0])) + ")")
        axarr[i].set_yticks([])
    # add rectangles denoting HA and NA on last subplot:
    # rectangles = []
    # axarr[-1].axhspan(ymin=0, ymax=10, xmin=0, xmax=0.55)
    # axarr[-1].axhspan(ymin=0, ymax=10, xmin=0.55, xmax=1.0, facecolor='r')
    # axarr[-1].set_xticks([])
    # axarr[-1].set_yticks([])
    # axarr[-1].axison = False
    filename = "".join([FIG_PATH, p,  "_responses_by_clusters_n_", str(num_clusters), ".png"])
    f.savefig(filename, dpi=200)    

    # add text on top of rectangles:
    # for r in arange(len(rectangles)):
    #     rx, ry = rectangles[r].get_xy()
    #     cx = rx + rectangles[r].get_width()/2.0
    #     cy = ry + rectangles[r].get_height()/2.0

    #     ax.annotate('BLA', (cx, cy), color='w', weight='bold', 
    #             fontsize=6, ha='center', va='center')
    # #rect = matplotlib.patches.Rectangle( (0, 0), width=clusters[p].len, height=10000)
    #axarr[-1].add_patch(rect)

# Plot median response per cluster
# for p in ['Vic_HA', 'Vic_NA']:
    
#     f, axarr = plt.subplots(num_clusters,1)    
#     f.set_tight_layout(True)
#     f.set_size_inches(18,11)
#     # plot clusters  
#     for i in np.arange(num_clusters):
        
#         axarr[i].bar(np.arange(len(ind_dict[p])), np.median(df[ind_dict[p]].loc[clusters[p] == i+1].T, axis=1))
#         axarr[i].set_title(p + " cluster " + str(i+1) + " (n = " + str(len(np.where([clusters[p] == i+1])[0])) + ")")
#         axarr[i].set_yticks([])
#         axarr[i].set_ylim(0, 30000)

#     filename = "".join([FIG_PATH, p,  "_median_responses_by_clusters_n_", str(num_clusters), ".png"])
#     f.savefig(filename, dpi=200)    

# Plot median response per pairs of clusters of interest
for p in ['Vic_HA']:
    
    ind = np.arange(len(ind_dict[p]))+1  # the x locations for the groups
    width = 0.3       # the width of the bars

    colors = ['r','b']
    # plot pairs of clusters with significant differences 
    for i in sig_df_HA.index:

        f, ax = plt.subplots()
        f.set_tight_layout(True)
        f.set_size_inches(18,11)

        #med_resp = np.median(df[ind_dict[p]].loc[clusters[p] == i[0]].T, axis=1)
        rects1 = ax.bar(ind, med_resp[:,i[0]-1], width, color=colors[0])
        
        #med_resp = np.median(df[ind_dict[p]].loc[clusters[p] == i[1]].T, axis=1)
        rects2 = ax.bar(ind+width, med_resp[:,i[1]-1], width, color=colors[1])
        ax.legend(i)
        ax.set_title(p + " median responses of clusters")
    
        filename = "".join([FIG_PATH, p,  "_bar_median_responses_clusters_", str(i[0]), "_",str(i[1]), ".png"])
        f.savefig(filename, dpi=200)    


# plot correlation matrix of HA and NA clustering on Victoria
# for p in ['Vic_HA', 'Vic_NA']:

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
for p in ['Vic_HA', 'Vic_NA']:
    colInd = hcp.plotHColCluster(df[ind_dict[p]].T,method='complete', metric='spearman',titleStr=p,vRange=(0.5, 1))
    f = gcf()
    filename = "".join([FIG_PATH + p + "_clustering_matrix_and_dendrogram"])
    f.savefig(filename, dpi=200)

# train classifier to predict infection status, or treatment assignment using scikit learn
Vic_data      = np.double(df[np.concatenate([ind_dict['Vic_HA'],ind_dict['Vic_NA']],axis=0)].as_matrix())
#treatment_labels = np.asarray(df['treatment_label']) 
#treatment_labels = np.asarray(df['case_control_label'])

HAI_values = df[Vic_assays[0]]
pos_label_inds = mlab.find(HAI_values > 40)
neg_label_inds = mlab.find(HAI_values < 40)
label_inds = np.concatenate((pos_label_inds,neg_label_inds),axis=0)

Vic_HAI_data = Vic_data[label_inds,:]
Vic_HAI_values = HAI_values[label_inds]
Vic_HAI_labels = np.zeros(shape=Vic_HAI_values.shape)
Vic_HAI_labels[mlab.find(Vic_HAI_values > 80)] = 1
Vic_HAI_labels[mlab.find(Vic_HAI_values < 40)] = 0


# here take quantile approach due to lack of known cutoff:
Vic_NT_values = df[Vic_assays[1]]
#NT_q_inds = pd.qcut(NT_values,q=4, labels=False)
# pos_label_inds = mlab.find(NT_q_inds ==0)
# neg_label_inds = mlab.find(NT_q_inds == 3)
pos_label_inds = mlab.find(Vic_NT_values > 1000)
neg_label_inds = mlab.find(Vic_NT_values <= 100)
label_inds = np.concatenate((pos_label_inds,neg_label_inds),axis=0)

Vic_NT_data   = Vic_data[label_inds,:]
Vic_NT_values = NT_values[label_inds]
Vic_NT_labels = np.zeros(shape=Vic_NT_values.shape)
Vic_NT_labels[mlab.find(Vic_NT_values > 1000)] = 1
Vic_NT_labels[mlab.find(Vic_NT_values <= 100)] = 0

for d in ['Vic_HAI', 'Vic_NT']:

    logistic = linear_model.LogisticRegression(C=0.1,penalty='l2')
    scores = []
    num_folds = 10
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    all_tpr = []

    if(d == 'Vic_HAI'):
        data   = Vic_HAI_data
        labels = Vic_HAI_labels
    elif(d == 'Vic_NT'):
        data   = Vic_NT_data
        labels = Vic_NT_labels


    k_fold = cross_validation.StratifiedKFold(labels, n_folds=num_folds)
    for k, (train, test) in enumerate(k_fold):
        logistic.fit(data[train], labels[train])
        curr_score = logistic.score(data[test], labels[test])
        scores.append(curr_score)
        print("[fold {0}], score: {1:.2f}".
          format(k, curr_score))

        probas_ = logistic.fit(data[train], labels[train]).predict_proba(data[test])
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = metrics.roc_curve(labels[test], probas_[:, 1])
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        roc_auc = metrics.auc(fpr, tpr)
        #plt.plot(fpr, tpr, lw=1, label='ROC fold %d (area = %0.2f)' % (k, roc_auc))

    f = figure();
    plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6))
    mean_tpr /= num_folds
    mean_tpr[-1] = 1.0
    mean_auc = metrics.auc(mean_fpr, mean_tpr)
    plt.plot(mean_fpr, mean_tpr, 'k--',lw=2)
    plt.title('Mean ROC (AUC = %0.2f)' % mean_auc)

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

    filename = "".join([FIG_PATH, d, "_ROC.png"])
    f.savefig(filename, dpi=200)   





# alphas = np.logspace(-4, -6, 30)
# num_folds = 10;
# k_fold = cross_validation.KFold(len(Vic_HA_data), n_folds=num_folds)
# best_alphas = np.zeros(shape=(num_folds))
# for k, (train, test) in enumerate(k_fold):
#     scores = []
#     for alpha in alphas:
#         lasso = linear_model.Lasso(alpha = alpha, max_iter=20, normalize=True)
#         lasso.fit(Vic_HA_data[train], treatment_labels[train])
#         pred_labels = lasso.predict(Vic_HA_data[train])
#         r, p = scipy.stats.pearsonr(pred_labels,treatment_labels[train])
#         scores.append(r)
#     max_ind = np.argmax(scores)
#     best_alphas[k]= alphas[max_ind]
#     print("[fold {0}],  alpha {1:.7f} r = {2:.2f},  p = {3:.2f}".
#         format(k, alphas[max_ind], scores[max_ind], p))

#     final_alpha = 1.00000000e-06
#     lasso = linear_model.Lasso(alpha = final_alpha, max_iter=100, normalize=True)
#     lasso.fit(Vic_HA_data[train], treatment_labels[train])
#     pred_labels = lasso.predict(Vic_HA_data[test])
#     r, p = scipy.stats.pearsonr(pred_labels,treatment_labels[test])
#     print("[fold {0}],  test error r = {1:.2f},  p = {2:.2f}".
#         format(k, r, p))





