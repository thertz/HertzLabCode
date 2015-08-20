from numpy import *
from pylab import *
from numpy.random import permutation,seed
import pandas as pd

__all__ = ['scatterdots',
           'myboxplot',
           'manyboxplots']

def scatterdots(data, x, width=0.8, returnx=False, rseed=820, **kwargs):
    """Dots plotted with random x-coordinates and y-coordinates from data array.

    Parameters
    ----------
    data : ndarray
    x : float
        Specifies the center of the dot cloud on the x-axis.
    width : float
        Specifies the range of the dots along the x-axis.
    returnx : bool
        If True, return the x-coordinates of the plotted data points.
    rseed : float
        Random seed. Defaults to a constant so that regenerated figures of
        the same data are identical.

    Returns
    -------
    Optionally returns the x-coordinates as plotted."""

    seed(rseed)

    if data is None or len(data) == 0:
        if returnx:
            return None
        return
    if not type(data) == ndarray:
        data = array(data)

    validi = arange(len(data))
    if any(isnan(data)):
        validi = where(logical_not(isnan(data)))[0]
    ploty = data[validi]

    if len(ploty) == 0:
        if returnx:
            return None
        return
    w = width
    plotx = permutation(linspace(-w/2.,w/2.,len(ploty))+x)
    scatter(plotx,ploty,**kwargs)
    
    if returnx:
        outx = nan*ones(data.shape)
        outx[validi] = plotx
        return outx

def myboxplot_by_labels(data, labels, positions=None, width=0.8, boxcolor='black', scatterwidth=0.6, dotcolor='red', altDotcolor='gray', **kwargs):
    """Make multiple boxplots from a single vector, using a label vector with scatterdots overlaid.

    Parameters
    ----------
    data : np.ndarray or pd.Series
    
    labels : int or strings
        label vector denoting the group to which each datapoint belongs
    
    positions : [floats | None (default)]
        Position of boxes along x-axis. if None, will position each box using its label value.
    
    width : float
        Width of the box.
    
    boxcolor : mpl color
    
    scatterwidth : float
        Width of the spread of the data points.
    
    dotcolor : mpl color
    
    altDotcolor : mpl color
        Specify the color of the data points that are not in the subset.
    
    """

    if type(data) is pd.Series:
        data = data.values

    # various formatting parameters passed on to boxplot function
    if not 's' in kwargs:
        kwargs['s'] = 20
    if not 'marker' in kwargs:
        kwargs['marker'] = 'o'
    if not 'linewidths' in kwargs:
        kwargs['linewidths'] = 0.5

    """Boxplot with dots overlaid""" 
    data_list = []
    uniq_l = np.unique(labels)
    for l in uniq_l:
        data_list.append(data[np.where(labels == l)[0]])


    bp = boxplot(data_list,widths=width,sym='')
    for element in bp.keys():
        for b in bp[element]:
            b.set_color(boxcolor)

    kwargs['c'] = dotcolor
    for i in arange(len(uniq_l)):
        subsetx = scatterdots(data[np.where(labels == uniq_l[i])[0]], x=i+1, width=scatterwidth, returnx=True, **kwargs)
    
    # plot correct labels on xticks:
    a = gca()
    a.set_xticklabels(uniq_l, rotation=45)
def myboxplot(data,x=1,width=0.8,boxcolor='black',scatterwidth=0.6,dotcolor='red',returnx=False,subsetInd=None,altDotcolor='gray',**kwargs):
    """Make a boxplot with scatterdots overlaid.

    Parameters
    ----------
    data : np.ndarray or pd.Series
    x : float
        Position of box along x-axis.
    width : float
        Width of the box.
    boxcolor : mpl color
    scatterwidth : float
        Width of the spread of the data points.
    dotcolor : mpl color
    subsetInd : boolean or int index
        Indicates a subset of the data that should be summarized in the boxplot.
        However, all data points will be plotted.
    altDotcolor : mpl color
        Specify the color of the data points that are not in the subset.
    returnx : bool
        Return the x-coordinates of the data points.

    Returns
    -------
    outx : np.ndarray
        Optionall, an array of the x-coordinates as plotted."""

    if type(data) is pd.Series:
        data = data.values

    if not subsetInd is None:
        if not (subsetInd.dtype == array([0,1],dtype=bool).dtype):
            tmp = zeros(data.shape,dtype=bool)
            tmp[subsetInd] = True
            subsetInd = tmp
    else:
        subsetInd = ones(data.shape,dtype=bool)

    if not 's' in kwargs:
        kwargs['s'] = 20
    if not 'marker' in kwargs:
        kwargs['marker'] = 'o'
    if not 'linewidths' in kwargs:
        kwargs['linewidths'] = 0.5

    """Boxplot with dots overlaid"""
    outx = zeros(data.shape)
    if subsetInd.sum() > 0:
        bp = boxplot(data[subsetInd],positions=[x],widths=width,sym='')
        for element in bp.keys():
            for b in bp[element]:
                b.set_color(boxcolor)

        kwargs['c'] = dotcolor
        subsetx = scatterdots(data[subsetInd], x=x, width=scatterwidth, returnx=True, **kwargs)
        outx[subsetInd] = subsetx

    if (~subsetInd).sum() > 0:
        kwargs['c'] = altDotcolor
        subsetx = scatterdots(data[~subsetInd], x=x, width=scatterwidth, returnx=True, **kwargs)
        outx[~subsetInd] = subsetx

    if returnx:
        return outx

def manyboxplots(df,cols=None,colLabels=None,annotation='N',horizontal=False,vRange=None,**kwargs):
    """Series of boxplots along x-axis (or flipped horizontally along y-axis [NOT IMPLEMENTED])

    Optionally add annotation for each boxplot with:
        (1) "N"
        (2) "pctpos" (response rate, by additionally specifying responders)
            NOT YET IMPLEMENTED

    Parameters
    ----------
    df : pd.DataFrame
    cols : list
        Column names to be plotted
    colLabels : list
        Column labels (optional)
    annotation : str or None
        Specifies what the annotation should be: "N" or "pctpos"
    horizontal : bool
        Specifies whether boxplots should be vertical (default, False) or horizontal (True)
    kwargs : additional arguments
        Passed to myboxplot function to specify colors etc."""
    if cols is None:
        cols = df.columns
    if colLabels is None:
        colLabels = cols
    elif len(colLabels)<cols:
        colLabels += cols[len(colLabels):]

    axh = gca()

    for x,c in enumerate(cols):
        myboxplot(df[c].dropna(), x=x, **kwargs)

    if not vRange is None:        
        ylim(vRange)
    yl = ylim()
    annotationKwargs = dict(xytext = (0,-10), textcoords = 'offset points', ha='center',va='top',size='medium')
    for x,c in enumerate(cols):
        tmp = df[c].dropna()

        if annotation == 'N':
            annotate('%d' % len(tmp),xy=(x,yl[1]), **annotationKwargs)
        elif annotation == 'pctpos':
            pass
        
    xlim((-1,x+1))
    xticks(arange(x+1))
    xlabelsL = axh.set_xticklabels(colLabels,fontsize='large',rotation=90,fontname='Consolas')
