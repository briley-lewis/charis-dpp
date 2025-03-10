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
from alive_progress import alive_bar
from skimage.measure import centroid

def register_cubes(datadir,offset=[31.758, 16.542],save_hypercube=False,save_plots=False,centroiding=True):
    """
    Function for registering CHARIS pol image cubes via cross-correlation. 
    Assumes data is in the format provided by CHARIS-DPP of n0001left.fits, n0001right.fits, etc.
    Saves registered data as "n0001left_reg.fits", "n0001right_reg.fits", etc. to be compatible with next steps in CHARIS-DPP.

    Parameters
    ----------
    datadir : string
        directory containing the data you want to register
        should have names like "n0001left.fits"
    offset : tuple, optional
        offset in pixels to shift the data by if centroiding=False, default is [31.758, 16.542] (from CHARIS-DPP charis_pdi_register_cube)
    save_hypercube : bool, optional
        whether to save the hypercube of registered data, default is False
    save_plots : bool, optional
        whether to save plots of the cross-correlation, default is False
    centroiding : bool, optional
        whether to use centroiding to determine the shift, default is True (False uses the offset parameter)
        
    Returns
    -------
    None
    """

    prep_files_right = glob.glob(datadir+'*right.fits')
    prep_files_left = glob.glob(datadir+'*left.fits')

    print(len(prep_files_left), 'files left', len(prep_files_right), 'files right')

    #create centered reference cubes for left and right -- this should be centroiding, not just a guess with the polcent_offset
    print("creating centered reference cubes")
    cntred_cubel = fits.getdata(prep_files_left[0])
    testcube_l = fits.getdata(prep_files_left[0])
    cent = np.shape(testcube_l)[1]/2
    for i in range(np.shape(testcube_l)[0]):
        #hardcoded values here are a hack to ignore junk at the edge of the FOV in centroiding
        mask = np.zeros(np.shape(testcube_l[i])) # create a mask with the image's shape
        mask[50:120,35:100] = 1 # create a mask
        offset_centroid_y = -centroid(testcube_l[i]*mask)[0]+cent
        offset_centroid_x = -centroid(testcube_l[i]*mask)[1]+cent
        cntred_image = shift(testcube_l[i], (offset_centroid_y,offset_centroid_x))
        cntred_cubel[i,:,:] = cntred_image

    cntred_cuber = fits.getdata(prep_files_right[0])
    testcube_r = fits.getdata(prep_files_right[0])
    for i in range(np.shape(testcube_r)[0]):
        #hardcoded values here are a hack to ignore junk at the edge of the FOV in centroiding
        mask = np.zeros(np.shape(testcube_r[i])) # create a mask with the image's shape
        mask[80:150,98:163] = 1 # create a mask 
        offset_centroid_x = -centroid(testcube_r[i]*mask)[1]+cent
        offset_centroid_y = -centroid(testcube_r[i]*mask)[0]+cent
        cntred_image = shift(testcube_r[i], (offset_centroid_y,offset_centroid_x))
        cntred_cuber[i,:,:] = cntred_image

    #create mask for aliasing effects
    mask_left = np.zeros(cntred_image.shape)
    polycorners_left = np.array([[(94,50),(145,50),(145,100),(115,158),(47,140)]])
    cv2.fillPoly(mask_left, polycorners_left,color=1)
    mask_left[mask_left==0] = np.nan

    mask_right = np.zeros(cntred_image.shape)
    polycorners_right = np.array([[(94,50),(145,50),(145,84),(115,146),(47,140)]])
    cv2.fillPoly(mask_right, polycorners_right,color=1)
    mask_right[mask_right==0] = np.nan

    #cross-correlation to center all cubes (left)
    print("starting cross-correlation for left cubes")
    count = 0
    cntred_data_hypercubel = np.zeros((len(prep_files_left),np.shape(testcube_l)[0],np.shape(testcube_l)[1],np.shape(testcube_l)[2]))
    with alive_bar(len(prep_files_left),title="Frames completed") as bar:
        for file in prep_files_left:
            ##read in cube
            data = fits.getdata(file)
            header = fits.getheader(file)
            ##create a copy to save centered data
            cntred_data_tmp = np.copy(data)
            for i in range(np.shape(testcube_l)[0]):
                ##calculate shift necessary to center
                shiftval, error, diffphase = phase_cross_correlation(
                cntred_cubel[i], data[i], upsample_factor=50)
                ##shift cube to reference center
                cntred_image = shift(data[i], shiftval)
                cntred_data_tmp[i,:,:] = cntred_image*mask_left
            header['HISTORY'] = 'Registered to data cube 0001left.fits for each wavelength slice'
            header['HISTORY'] = 'Registration completed with skimage.phase_cross_correlation bfrom nonstandard_pdi_reg by B. Lewis (2025)'
            fits.writeto(datadir+'n{:04d}left_reg.fits'.format(count+1),cntred_data_tmp,header,overwrite=True)
            cntred_data_hypercubel[count,:,:,:] = cntred_data_tmp
            count+=1
            bar()

    print("left cross-correlation complete")

    #cross-correlation to center all cubes (right)
    print("starting cross-correlation for right cubes")
    count = 0
    cntred_data_hypercuber = np.zeros((len(prep_files_right),np.shape(testcube_l)[0],np.shape(testcube_l)[1],np.shape(testcube_l)[2]))
    with alive_bar(len(prep_files_right),title="Frames completed") as bar:
        for file in prep_files_right:
            ##read in cube
            data = fits.getdata(file)
            header = fits.getheader(file)
            ##create a copy to save centered data
            cntred_data_tmp = np.copy(data)
            for i in range(np.shape(testcube_r)[0]):
                ##calculate shift necessary to center
                shiftval, error, diffphase = phase_cross_correlation(
                cntred_cuber[i], data[i], upsample_factor=50)
                ##shift cube to reference center
                cntred_image = shift(data[i], shiftval)
                cntred_data_tmp[i,:,:] = cntred_image*mask_right
            header['HISTORY'] = 'Registered to data cube 0001right.fits for each wavelength slice'
            header['HISTORY'] = 'Registration completed with skimage.phase_cross_correlation bfrom nonstandard_pdi_reg by B. Lewis (2025)'
            fits.writeto(datadir+'n{:04d}right_reg.fits'.format(count+1),cntred_data_tmp,header,overwrite=True)
            cntred_data_hypercuber[count,:,:,:] = cntred_data_tmp
            count+=1
            bar()

    print("right cross-correlation complete")

    #plot left crosscorr
    fig = plt.figure(figsize=(8, 3))
    ax1 = plt.subplot(1, 3, 1)
    ax2 = plt.subplot(1, 3, 2, sharex=ax1, sharey=ax1)
    ax3 = plt.subplot(1, 3, 3)

    ax1.imshow(cntred_cubel[19,:,:], cmap='gray',vmax=1e3)
    ax1.set_axis_off()
    #ax1.scatter(100,100,marker='x',color='red')
    ax1.set_title('Reference image')

    ax2.imshow(testcube_l[0], cmap='gray',vmax=1e3)
    ax2.set_axis_off()
    ax2.set_title('Offset image')

    ax3.imshow(cntred_data_hypercubel[8,19,:,:], cmap='gray',vmax=1e3)
    ax3.set_axis_off()
    ax3.set_title('Centered data image')
    #ax3.scatter(100,100,marker='x',color='red')

    plt.suptitle('Left Pol Image Cross Correlation')
    if save_plots:
        plt.savefig(datadir+'left-crosscorr.png')
    plt.show()

    ##plot right cross corr
    fig = plt.figure(figsize=(8, 3))
    ax1 = plt.subplot(1, 3, 1)
    ax2 = plt.subplot(1, 3, 2, sharex=ax1, sharey=ax1)
    ax3 = plt.subplot(1, 3, 3)

    ax1.imshow(cntred_cuber[0,:,:], cmap='gray',vmax=10e3)
    ax1.set_axis_off()
    ax1.scatter(100,100,marker='x',color='red')
    ax1.set_title('Reference image')

    ax2.imshow(testcube_r[0], cmap='gray',vmax=10e3)
    ax2.set_axis_off()
    ax2.set_title('Offset image')

    ax3.imshow(cntred_data_hypercuber[0,19,:,:], cmap='gray',vmax=1e3)
    #ax3.scatter(100,100,marker='x',color='red')
    ax3.set_axis_off()
    ax3.set_title('Centered data image')

    plt.suptitle('Right Pol Image Cross Correlation')
    if save_plots:
        plt.savefig(datadir+'right-crosscorr.png')
    plt.show()

    if save_hypercube:
        hdu = fits.PrimaryHDU(data=cntred_data_hypercubel)
        hdu.writeto(datadir+'centered_hypercube_left.fits',overwrite=True)
        hdu = fits.PrimaryHDU(data=cntred_data_hypercuber)
        hdu.writeto(datadir+'centered_hypercube_right.fits',overwrite=True)

    print('Image registration complete for',len(prep_files_left),'frames')

    return None

import sys
if __name__ == '__main__':
    if len(sys.argv)<2:
        datadir = input('Please specify a directory (full path) containing the data you want to register \n')
    else:
        datadir = sys.argv[1]
    offset_new = input('Please specify the offset you want to use for registration in the form [x,y], or press enter to use the default [31.758, 16.542] \n')
    if offset_new != '':
        offset = offset_new
        register_cubes(datadir,offset=offset_new)
    else:
        register_cubes(datadir)