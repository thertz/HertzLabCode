""" utiltlity functions for antigen array data Hertz Lab"""
import pandas as pd
import matplotlib.mlab as mlab
import scipy.stats
import numpy as np
import itertools
import ipdb as ipdb

__all__ = ['']

def compare_response_peaks(resp_df, peak_threshold=30000):
    """
    compare the responses of peaks in the resp_df pandas dataframe (antigen array data).

    Parameters:
    ----------
    resp_df: pandas.DataFrame
        dataframe of responses of the samples to compare
    peak_threshold: int
        threshold defining response peaks. Default set to 30,000.
    
    Returns:
    -------
        output of responses of all samples to peaks in each sample in the dataframe
    """
    resp_df = 