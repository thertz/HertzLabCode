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
SAVE_PATH = ROOT_PATH + 'Data/DeImm/'
FIG_PATH = SAVE_PATH + 'Figs/'

prot_paths = {'Fabri': 'Fabri', 'Factor7': 'Factor7/WT', 'Factor7_Mut': 'Factor7/MUT/'}

percentile_threshold = 2 # threshold by which to define an epitope

allele_frequencies = {'HLA-DRB1*01:01': 0.054, 'HLA-DRB1*03:01': 0.137, 'HLA-DRB1*04:01': 0.046,
                      'HLA-DRB1*04:05': 0.062, 'HLA-DRB1*07:01': 0.135, 'HLA-DRB1*08:02': 0.049, 
                      'HLA-DRB1*09:01': 0.062, 'HLA-DRB1*11:01': 0.118, 'HLA-DRB1*12:01': 0.039,  
                      'HLA-DRB1*13:02': 0.077, 'HLA-DRB1*15:01': 0.122, 'HLA-DRB3*01:01': 0.261,
                      'HLA-DRB3*02:02': 0.343, 'HLA-DRB4*01:01': 0.418, 'HLA-DRB5*01:01': 0.16}

IEDB_dict = {}
epitope_maps = {}
pop_epitope_map = {}
pop_epitope_map_weighted = {}

for prot in prot_paths.keys():
    
    IEDB_dict[prot] = {}
    epitope_maps[prot] = {}


    for fni, fileName in enumerate(glob.glob(SAVE_PATH + prot_paths[prot] + '/*.csv')):
        print(fni, fileName)
        curr_df = pd.read_csv(fileName)
        curr_df = curr_df.sort(columns='start')
        curr_df['is_epitope'] = curr_df.percentile_rank.map(lambda p: 1 if p<=percentile_threshold else 0)

        curr_allele = curr_df.iloc[0].allele
        
        IEDB_dict[prot][curr_allele] = curr_df
        epitope_maps[prot][curr_allele] = curr_df.is_epitope

        # init population epitope map only on first allele:
        if fni == 0:
            pop_epitope_map[prot] = np.zeros(len(epitope_maps[prot][curr_allele]))
            pop_epitope_map_weighted[prot] = np.zeros(len(epitope_maps[prot][curr_allele]))

        pop_epitope_map[prot] += epitope_maps[prot][curr_allele]
        pop_epitope_map_weighted[prot] += epitope_maps[prot][curr_allele]*allele_frequencies[curr_allele]

    # #plot population based epitope map:
    # f, axarr = plt.subplots(2,1)
    # f.set_tight_layout(True)
    # f.set_size_inches(18,11)
    # axarr[0].plot(pop_epitope_map[prot])
    # axarr[0].set_title("".join([prot, ' population based epitope map']))
    # axarr[0].set_xlabel('Epitope Start Position')
    # axarr[0].set_ylabel('# of epitopes')

    # axarr[1].plot(pop_epitope_map_weighted[prot]/2)
    # axarr[1].set_title("".join([prot, ' weighted population based epitope map']))
    # axarr[1].set_xlabel('Epitope Start Position')
    # axarr[1].set_ylabel('population coverage (%)')

    # filename = "".join([FIG_PATH, prot, "_epitope_maps_for_HLA_class_II_th_", str(percentile_threshold), ".png"])
    # f.savefig(filename, dpi=200)



    f, axarr = plt.subplots(1,1)
    f.set_tight_layout(True)
    f.set_size_inches(18,11)

    axarr.plot(pop_epitope_map_weighted[prot]/2)
    axarr.set_title("".join([prot, ' weighted population based epitope map']))
    axarr.set_xlabel('Epitope Start Position')
    axarr.set_ylabel('population coverage (%)')

    filename = "".join([FIG_PATH, prot, "Weighted_epitope_maps_for_HLA_class_II_th_", str(percentile_threshold), ".png"])
    f.savefig(filename, dpi=200)
