"""
Initialize array parameters All of these can be modified as required.
"""
#-------------------------------------------------------------------------------------------------------------------------------------
# project specific parameters - this is for H7N9, not relevant for other projects
#-------------------------------------------------------------------------------------------------------------------------------------




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
bg_threshold = 2000 # threshold by which to flag antigens as ones that have high background and should be removed from analysis (typically using a BSA negative control).
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


