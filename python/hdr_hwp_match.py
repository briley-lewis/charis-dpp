"""HWP matching using information in CHARIS headers. Generally, use the CHARIS-DPP, 
but in case that script doesn't work / produces nonsense, here's another thing to try.

Written by B. Lewis (2025)"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob
from alive_progress import alive_bar

def hwp_match(datadir):
        """
    Function for getting HWP angle positions from CHARIS headers

    Parameters
    ----------
    datadir : string
        directory containing the "prep" folder of data you want to register
        e.g. /Users/brileylewis/Documents/postdoc/europa/2023-07-30/reduc/
        should have names like "n0001left.fits"
        
    Returns
    -------
    None"""
    print('work in progress, no hwp log for you yet')

import sys
if __name__ == '__main__':
    if len(sys.argv)<2:
        datadir = input('Please specify a directory (full path) containing the data for which you want to get HWP angles \n')
    else:
        datadir = sys.argv[1]
    hwp_match(datadir)