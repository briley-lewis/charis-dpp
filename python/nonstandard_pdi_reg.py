"""Image registration for CHARIS pol data, where data is in a non-standard format (e.g. no satellite spots).

Requires data structured in folders as in the CHARIS pipeline, returns data such that you can continue to next
steps with CHARIS-DPP seamlessly. Must run charis_pdi_split_pols.pro before this script.

Written by B. Lewis (2025)"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob
from skimage import data
from skimage.registration import phase_cross_correlation
from skimage.registration._phase_cross_correlation import _upsampled_dft
from scipy.ndimage import fourier_shift, shift
import cv2

def register_cubes(datadir,offset=[0,0]):
    """docstring go here"""

    prep_files_right = glob.glob(datadir+'*right.fits')
    prep_files_left = glob.glob(datadir'*left.fits')

    print(len(prep_files_left), 'files left', len(prep_files_right), 'files right')

    # Nominal center values taken from CHARIS-DPP charis_pdi_register_cube
    xcent_right_nom = 0
    ycent_right_nom = 0
    xcent_left_nom = 0
    ycent_left_nom = 0

    ###something go here to check if offset center looks correct