import palettable
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec
from matplotlib import cm
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch
import numpy as np
import matplotlib as mpl
import pylab
import itertools



__all__ = ['plotHCluster','plotHColCluster','plotCorrHeatmap','mapColors2Labels']

def clean_axis(ax):
    """Remove ticks, tick labels, and frame from axis"""
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

def mapColors2Labels(labels,setStr='Set1'):
    """Return pd.Series of colors based on labels"""
    N = max(3,min(12,len(np.unique(labels))))
    cmap_name = "palettable.colorbrewer.qualitative." + setStr + "_" + str(N) + ".mpl_colors"
    cmap=eval(cmap_name)
    cmapLookup = {k:col for k,col in zip(sorted(np.unique(labels)),itertools.cycle(cmap))}
    return labels.map(cmapLookup.get)

def computeDMat(df,metric,minN=1):
    if metric=='spearman':
        dmat = 1 - df.T.corr(method='spearman',min_periods = minN).values
        dmat[np.isnan(dmat)]=0
    else:
        dmat = distance.squareform(distance.pdist(df,metric=metric))
    return dmat

def computeHCluster(df,method='complete',metric='euclidean',dmat=None,minN=1):
    """Compute dmat, clusters and dendrogram of df using
    the linkage method and distance metric given"""
    if dmat is None:
        dmat = computeDMat(df,metric,minN=minN)
    clusters = sch.linkage(dmat, method=method)
    den = sch.dendrogram(clusters,color_threshold=np.inf,no_plot=True)

    if metric in ['spearman','pearson','correlation']:
        """Translate correlation-based distances back to [-1,1] scale"""
        dmat = 1 - dmat
    return dmat, clusters, den

def testData(rows=50,columns=20):
    data=np.random.multivariate_normal(rand(columns),rand(columns,columns),rows)
    df=pd.DataFrame(data,columns=[''.join([lett]*9) for lett in 'ABCDEFGHIJKLMNOPQRST'])
    rowLabels=pd.Series(rand(rows).round(),index=df.index)
    columnLabels=pd.Series(rand(columns).round(),index=df.columns)
    return {'df':df,'row_labels':rowLabels,'col_labels':columnLabels}

def addColorbar(fig,cb_ax,data_ax,label='Correlation'):
    """Colorbar"""
    cb = fig.colorbar(data_ax,cb_ax) # note that we could pass the norm explicitly with norm=my_norm
    cb.set_label(label)
    """Make colorbar labels smaller"""
    for t in cb.ax.yaxis.get_ticklabels():
        t.set_fontsize('small')

def plotCorrHeatmap(df=None,metric='spearman',colInd=None,col_labels=None,titleStr=None,vRange=None,tickSz='small',cmap=None,dmat=None,cbLabel='Correlation',minN=1):
    """Plot a heatmap of a column-wise distance matrix defined by metric (can be 'spearman' as well)
    Can provide dmat as a pd.DataFrame instead of df.
    Optionally supply a column index colInd to reorder the columns to match a previous clustering
    Optionally, col_labels will define a color strip along the yaxis to show groups"""

    fig = plt.gcf()
    fig.clf()

    if df is None:
        """If df is None then data will come from a DataFrame dmat. Create dummy df for labels though"""
        df = pd.DataFrame(data=np.zeros((5,dmat.shape[0])),columns=dmat.columns)
        dmat = dmat.values

    if dmat is None:
        dmat = computeDMat(df.T,metric,minN=minN)

    if cmap is None:
        cmap=cm.RdBu_r

    if colInd is None:
        colInd=np.arange(dmat.shape[1])

    if col_labels is None:
        heatmapAX = fig.add_subplot(GridSpec(1,1,left=0.05,bottom=0.05,right=0.78,top=0.85)[0,0])
        scale_cbAX = fig.add_subplot(GridSpec(1,1,left=0.87,bottom=0.05,right=0.93,top=0.85)[0,0])
    else:
        col_cbAX = fig.add_subplot(GridSpec(1,1,left=0.05,bottom=0.05,right=0.08,top=0.85)[0,0])
        heatmapAX = fig.add_subplot(GridSpec(1,1,left=0.11,bottom=0.05,right=0.78,top=0.85)[0,0])
        scale_cbAX = fig.add_subplot(GridSpec(1,1,left=0.87,bottom=0.05,right=0.93,top=0.85)[0,0])
    
    if vRange is None:
        vmin,vmax = (-1,1)
        #vmin = dmat.flatten().min()
        #vmax = dmat.flatten().max()
    else:
        vmin,vmax=vRange
    my_norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    """Column label colorbar but along the rows"""
    if not col_labels is None:
        col_cbSE = mapColors2Labels(col_labels)
        col_axi = col_cbAX.imshow([[x] for x in col_cbSE.iloc[colInd].values],interpolation='nearest',aspect='auto',origin='lower')
        clean_axis(col_cbAX)

    """Heatmap plot"""
    axi = heatmapAX.imshow(dmat[colInd,:][:,colInd],interpolation='nearest',aspect='auto',origin='lower',norm=my_norm,cmap=cmap)
    clean_axis(heatmapAX)

    """Column tick labels along the rows"""
    if tickSz is None:
        heatmapAX.set_yticks([])
        heatmapAX.set_xticks([])
    else:
        heatmapAX.set_yticks(np.arange(dmat.shape[1]))
        heatmapAX.yaxis.set_ticks_position('right')
        heatmapAX.set_yticklabels(df.columns[colInd],fontsize=tickSz,fontname='Consolas')

        """Column tick labels"""
        heatmapAX.set_xticks(np.arange(dmat.shape[1]))
        heatmapAX.xaxis.set_ticks_position('top')
        xlabelsL = heatmapAX.set_xticklabels(df.columns[colInd],fontsize=tickSz,rotation=90,fontname='Consolas')

        """Remove the tick lines"""
        for l in heatmapAX.get_xticklines() + heatmapAX.get_yticklines(): 
            l.set_markersize(0)

    addColorbar(fig,scale_cbAX,axi,label=cbLabel)
    
    """Add title as xaxis label"""
    if not titleStr is None:
        heatmapAX.set_xlabel(titleStr,size='x-large')

def plotHColCluster(df,method='complete', metric='euclidean', col_labels=None,titleStr=None,vRange=None,tickSz='small',cmap=None,col_dmat=None,minN=1):
    """Perform hierarchical clustering on df columns and plot square heatmap of pairwise distances"""
    if cmap is None:
        if metric in ['spearman','pearson','correlation']:
            cmap=cm.RdBu_r
        else:
            cmap=cm.YlOrRd
    if metric in ['spearman','pearson','correlation']:
        colorbarLabel = 'Correlation coefficient'
    else:
        colorbarLabel = ''


    fig = plt.gcf()
    fig.clf()
    #heatmapGS = gridspec.GridSpec(1,4,wspace=0.0,width_ratios=[0.25,0.01,2,0.15])
    if col_labels is None:
        """heatmapGS = gridspec.GridSpec(2,3,width_ratios=[0.2,2,0.15],height_ratios=[0.1,1])
        col_denAX = fig.add_subplot(heatmapGS[1,0])
        heatmapAX = fig.add_subplot(heatmapGS[1,1])
        scale_cbAX = fig.add_subplot(heatmapGS[1,2])"""
        col_denAX = fig.add_subplot(GridSpec(1,1,left=0.05,bottom=0.05,right=0.15,top=0.85)[0,0])
        heatmapAX = fig.add_subplot(GridSpec(1,1,left=0.16,bottom=0.05,right=0.78,top=0.85)[0,0])
        scale_cbAX = fig.add_subplot(GridSpec(1,1,left=0.87,bottom=0.05,right=0.93,top=0.85)[0,0])
    else:
        """heatmapGS = gridspec.GridSpec(2,3,width_ratios=[0.2,0.01,2,0.15],height_ratios=[0.1,1])
        col_denAX = fig.add_subplot(heatmapGS[1,0])
        col_cbAX = fig.add_subplot(heatmapGS[1,1])
        heatmapAX = fig.add_subplot(heatmapGS[1,2])
        scale_cbAX = fig.add_subplot(heatmapGS[1,3])"""

        """TODO: work on row_cbAX so that I can have the data labels on the top and left"""
        col_denAX = fig.add_subplot(GridSpec(1,1,left=0.05,bottom=0.05,right=0.15,top=0.85)[0,0])
        col_cbAX = fig.add_subplot(GridSpec(1,1,left=0.16,bottom=0.05,right=0.19,top=0.85)[0,0])
        #row_cbAX = fig.add_subplot(GridSpec(1,1,left=0.2,bottom=0.83,right=0.78,top=0.86)[0,0])
        heatmapAX = fig.add_subplot(GridSpec(1,1,left=0.2,bottom=0.05,right=0.78,top=0.85)[0,0])
        scale_cbAX = fig.add_subplot(GridSpec(1,1,left=0.87,bottom=0.05,right=0.93,top=0.85)[0,0])
    
    col_dmat,col_clusters,col_den=computeHCluster(df.T,method,metric,col_dmat,minN=minN)

    if vRange is None:
        if metric in ['spearman','pearson','correlation']:
            vmin,vmax = (-1,1)
        else:
            vmin = col_dmat.flatten().min()
            vmax = col_dmat.flatten().max()
    else:
        vmin,vmax=vRange
    my_norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    """Column dendrogaram but along the rows"""
    pylab.axes(col_denAX)
    col_denD = sch.dendrogram(col_clusters,color_threshold=np.inf,orientation='right', count_sort=True, labels=col_labels)
    colInd = col_denD['leaves']
    clean_axis(col_denAX)

    """Column label colorbar but along the rows"""
    if not col_labels is None:
        col_cbSE = mapColors2Labels(col_labels)
        col_axi = col_cbAX.imshow([[x] for x in col_cbSE.iloc[palettable].values],interpolation='nearest',aspect='auto',origin='lower')
        clean_axis(col_cbAX)
        #row_axi = row_cbAX.imshow([x for x in col_cbSE.iloc[colInd].values],interpolation='nearest',aspect='auto',origin='lower')
        #clean_axis(row_cbAX)

    """Heatmap plot"""
    axi = heatmapAX.imshow(col_dmat[colInd,:][:,colInd],interpolation='nearest',aspect='auto',origin='lower',norm=my_norm,cmap=cmap)
    clean_axis(heatmapAX)

    """Column tick labels along the rows"""
    if tickSz is None:
        heatmapAX.set_yticks(())
        heatmapAX.set_xticks(())
    else:
        heatmapAX.set_yticks(np.arange(df.shape[1]))
        heatmapAX.yaxis.set_ticks_position('right')
        heatmapAX.set_yticklabels(df.columns[colInd],fontsize=tickSz,fontname='Consolas')

        """Column tick labels"""
        heatmapAX.set_xticks(np.arange(df.shape[1]))
        heatmapAX.xaxis.set_ticks_position('top')
        xlabelsL = heatmapAX.set_xticklabels(df.columns[colInd],fontsize=tickSz,rotation=90,fontname='Consolas')

        """Remove the tick lines"""
        for l in heatmapAX.get_xticklines() + heatmapAX.get_yticklines(): 
            l.set_markersize(0)

    addColorbar(fig,scale_cbAX,axi,label=colorbarLabel)

    """Add title as xaxis label"""
    if not titleStr is None:
        heatmapAX.set_xlabel(titleStr,size='x-large')
    return colInd

def plotHCluster(df, method='complete', metric='euclidean', clusterBool=[True,True],row_labels=None, col_labels=None, vRange=None,titleStr=None,xTickSz='small',yTickSz='small',cmap=None,minN=1):
    """Perform hierarchical clustering on df data columns (and rows) and plot results as
    dendrograms and heatmap.

    df - pd.DataFrame(), will use index and column labels as tick labels
    method and metric - parameters passed to scipy.spatial.distance.pdist and scipy.cluster.hierarchy.linkage
    row_labels - pd.Series with index same as df with values indicating groups (optional)
    col_labels - pd.Series with index same as columns in df with values indicating groups (optional)
    vMinMax - optional scaling, [vmin, vmax] can be derived from data
    clusterBool - [row, col] bool indicating whether to cluster along that axis
    """
    if cmap is None:
        cmap = cm.RdBu_r

    if vRange is None:
        vmin = df.min().min()
        vmax = df.max().max()
    else:
        vmin,vmax = vRange
    my_norm = mpl.colors.Normalize(vmin, vmax)

    fig = plt.gcf()
    fig.clf()
    heatmapGS = gridspec.GridSpec(3,3,wspace=0.0,hspace=0.0,width_ratios=[0.15,0.02,1],height_ratios=[0.15,0.02,1])

    if clusterBool[0]:
        row_dmat,row_clusters,row_den=computeHCluster(df,method,metric,minN=minN)

        """Dendrogarams"""
        row_denAX = fig.add_subplot(heatmapGS[2,0])
        row_denD = sch.dendrogram(row_clusters,color_threshold=np.inf,orientation='right')
        clean_axis(row_denAX)

        """Row colorbar"""
        if not row_labels is None:
            """NOTE: row_labels will not be index aware and must be in identical order as data"""
            row_cbSE = mapColors2Labels(row_labels,'Set1')
            row_cbAX = fig.add_subplot(heatmapGS[2,1])

            row_axi = row_cbAX.imshow([[x] for x in row_cbSE.iloc[row_denD['leaves']].values ],interpolation='nearest',aspect='auto',origin='lower')
            clean_axis(row_cbAX)
        rowInd=row_denD['leaves']
    else:
        rowInd=np.arange(df.shape[0])

    if clusterBool[1]:
        col_dmat,col_clusters,col_den=computeHCluster(df.T,method,metric,minN=minN)

        """Dendrogarams"""
        col_denAX = fig.add_subplot(heatmapGS[0,2])
        col_denD = sch.dendrogram(col_clusters,color_threshold=np.inf)
        clean_axis(col_denAX)

        """Column colorbar"""
        if not col_labels is None:
            col_cbSE = mapColors2Labels(col_labels)
            col_cbAX = fig.add_subplot(heatmapGS[1,2])
            col_axi = col_cbAX.imshow([list(col_cbSE.iloc[col_denD['leaves']])],interpolation='nearest',aspect='auto',origin='lower')
            clean_axis(col_cbAX)
        colInd=col_denD['leaves']
    else:
        colInd=np.arange(df.shape[1])
    
    """Heatmap plot"""
    heatmapAX = fig.add_subplot(heatmapGS[2,2])
    axi = heatmapAX.imshow(df.iloc[rowInd,colInd],interpolation='nearest',aspect='auto',origin='lower',norm=my_norm,cmap=cmap)
    clean_axis(heatmapAX)
    heatmapAX.grid(True)

    """Row tick labels"""
    heatmapAX.set_yticks(np.arange(df.shape[0]))
    ylabelsL = None
    if not yTickSz is None:
        heatmapAX.yaxis.set_ticks_position('right')
        ylabelsL = heatmapAX.set_yticklabels(df.index[rowInd],fontsize=yTickSz,fontname='Consolas')
    else:
        ylabelsL = heatmapAX.set_yticklabels([])

    """Add title as xaxis label"""
    if not titleStr is None:
        heatmapAX.set_xlabel(titleStr,size='x-large')

    """Column tick labels"""
    heatmapAX.set_xticks(np.arange(df.shape[1]))
    xlabelsL = None
    if not xTickSz is None:
        xlabelsL = heatmapAX.set_xticklabels(df.columns[colInd],fontsize=xTickSz,rotation=90,fontname='Consolas')

    """Remove the tick lines"""
    for l in heatmapAX.get_xticklines() + heatmapAX.get_yticklines(): 
        l.set_markersize(0)

    """Colorbar"""
    scaleGS = gridspec.GridSpec(10,15,wspace=0.,hspace=0.)
    scale_cbAX = fig.add_subplot(scaleGS[:2,0]) # colorbar for scale in upper left corner
    cb = fig.colorbar(axi,scale_cbAX) # note that we could pass the norm explicitly with norm=my_norm
    cb.set_label('Measurements')
    cb.ax.yaxis.set_ticks_position('left') # move ticks to left side of colorbar to avoid problems with tight_layout
    cb.ax.yaxis.set_label_position('left') # move label to left side of colorbar to avoid problems with tight_layout
    #cb.outline.set_linewidth(0)
    """Make colorbar labels smaller"""
    for t in cb.ax.yaxis.get_ticklabels():
        t.set_fontsize('small')
    scaleGS.tight_layout(fig,h_pad=0.0,w_pad=0.0)

    heatmapGS.tight_layout(fig,h_pad=0.1,w_pad=0.5)

    handles = dict(cb=cb, heatmapAX=heatmapAX, fig=fig,xlabelsL=xlabelsL,ylabelsL=ylabelsL,heatmapGS=heatmapGS)
    return rowInd,colInd,handles

def plotDfHeatmap(df,vRange=None,cmap=None,tickSz='large',titleStr=None,cbLabel=None):
    """Simple heatmap plot of data in df"""
    fig = plt.gcf()
    fig.clf()

    heatmapAX = fig.add_subplot(GridSpec(1,1,left=0.11,bottom=0.05,right=0.78,top=0.85)[0,0])
    scale_cbAX = fig.add_subplot(GridSpec(1,1,left=0.87,bottom=0.05,right=0.93,top=0.85)[0,0])
    
    if cmap is None:
        cmap=cm.RdBu_r
    if cbLabel is None:
        cbLabel=''
    if vRange is None:
        vmin = df.values.flatten().min()
        vmax = df.values.flatten().max()
    else:
        vmin,vmax=vRange
    my_norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    """Heatmap plot"""
    axi = heatmapAX.imshow(df.values.astype(float),interpolation='nearest',aspect='auto',origin='lower',norm=my_norm,cmap=cmap)
    clean_axis(heatmapAX)

    """Column tick labels along the rows"""
    if tickSz is None:
        heatmapAX.set_yticks([])
        heatmapAX.set_xticks([])
    else:
        heatmapAX.set_yticks(np.arange(df.shape[0]))
        heatmapAX.yaxis.set_ticks_position('right')
        heatmapAX.set_yticklabels(df.index,fontsize=tickSz,fontname='Consolas')

        """Column tick labels"""
        heatmapAX.set_xticks(np.arange(df.shape[1]))
        heatmapAX.xaxis.set_ticks_position('top')
        xlabelsL = heatmapAX.set_xticklabels(df.columns,fontsize=tickSz,rotation=90,fontname='Consolas')

        """Remove the tick lines"""
        for l in heatmapAX.get_xticklines() + heatmapAX.get_yticklines(): 
            l.set_markersize(0)

    addColorbar(fig,scale_cbAX,axi,label=cbLabel)
    
    """Add title as xaxis label"""
    if not titleStr is None:
        heatmapAX.set_xlabel(titleStr,size='x-large')

def normalizeAxis(df,axis=0,useMedian=False):
    """Normalize along the specified axis by
    subtracting the mean and dividing by the stdev.

    Uses df functions that ignore NAs

    Parameters
    ----------
    df : pd.DataFrame
    axis : int
        Normalization along this axis. (e.g. df.mean(axis=axis))

    Returns
    -------
    out : pd.DataFrame"""

    tmp = df.copy()
    retile = ones(len(df.shape))
    retile[axis] = df.shape[axis]
    if useMedian:
        tmp = tmp - tile(tmp.median(axis=axis).values,retile)
    else:
        tmp = tmp - tile(tmp.mean(axis=axis).values,retile)
    tmp = tmp / tile(tmp.std(axis=axis).values,retile)
    return tmp
