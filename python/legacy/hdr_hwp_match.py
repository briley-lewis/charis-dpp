"""HWP matching using information in CHARIS headers. Generally, use the CHARIS-DPP, 
but in case that script doesn't work / produces nonsense, here's another thing to try.

Written by B. Lewis (2025)"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob
#from alive_progress import alive_bar

def hwp_match(f_sat):
    """
    Function for getting HWP angle positions from CHARIS headers.
    Not as flexible as the pipeline, but a barebones option for if things get weird

    Parameters
    ----------
    f_sat: string
        numbers of files to parse for hwp matching - should be a list of two numbers, the first and last files to parse written as 5,10
        
    Returns
    -------
    None"""
    
    reducdir='./reduc/'
    datadir=reducdir+'reg/'
    hwpfileout='hwp_info.txt'

    f_sat = [int(f_sat.split(',')[0]),int(f_sat.split(',')[1])]

    # Get list of files
    files = glob.glob(datadir+'*reg.fits')
    if len(files)==0:
        print('No files found in directory '+datadir)
        print('trying cal suffix')
        files = glob.glob(datadir+'*cal.fits')
    if len(files)==0:
        raise ValueError('No files found in directory '+datadir+' Check your paths and try again')

    f_raw = sorted(glob.glob('./data/raw/*cube.fits'))

    hwp_ang = []
    hwp_pos = []
    hwp_cycles = []
    hwp_files = []

    #get file number from filename
    #get header angle position
    #set position number -- 0, 2, 1, 3 = 0, 45, 22.5, 67.5
    #set cycle number?? every 4, increments up at 0
    #exception for hwp ang = -1

    cycle_counter = -1
    for f_num in range(f_sat[0],f_sat[1]+1):
        #get hwp angle and set position
        hdul = fits.open(f_raw[f_num-1])
        #print(f_num)
        #print(f_raw[f_num])
        try:
            hwp_ang_i = hdul[3].header['RET-ANG1']
        except:
            hwp_ang_i = -1
            print('frame '+str(f_num)+' has no HWP angle')
        #print(hwp_ang_i)
        if hwp_ang_i==0:
            hwp_pos_i = 0
            cycle_counter += 1
        elif hwp_ang_i==45:
            hwp_pos_i = 2
        elif hwp_ang_i==22.5:
            hwp_pos_i = 1
        elif hwp_ang_i==67.5:
            hwp_pos_i = 3
        
        if hwp_cycles.count(cycle_counter)>=4:
            cycle_counter += 1
        
        if hwp_ang_i==-1:
            pass
        else:
            hwp_ang.append(hwp_ang_i)
            hwp_pos.append(hwp_pos_i)
            hwp_files.append(f_num)
            hwp_cycles.append(cycle_counter)
        hdul.close()

    #now, do a check each cycle index appears exactly four times
    cycles_to_remove = []
    for v in np.unique(hwp_cycles):
        if hwp_cycles.count(v)!=4:
            print('Cycle '+str(v)+' does not have 4 entries')
            cycles_to_remove.append(v)

    hwp_ang = np.asarray(hwp_ang)  
    hwp_ang = np.delete(hwp_ang, np.where(np.isin(hwp_cycles,cycles_to_remove)))
    hwp_pos = np.asarray(hwp_pos)
    hwp_pos = np.delete(hwp_pos, np.where(np.isin(hwp_cycles,cycles_to_remove)))
    hwp_files = np.asarray(hwp_files)
    hwp_files = np.delete(hwp_files, np.where(np.isin(hwp_cycles,cycles_to_remove)))
    hwp_cycles = np.asarray(hwp_cycles)
    hwp_cycles = np.delete(hwp_cycles, np.where(np.isin(hwp_cycles,cycles_to_remove)))

    print('incomplete cycles removed')

    hwp_info = np.column_stack((hwp_files, hwp_pos, hwp_ang, hwp_cycles))
    np.savetxt(hwpfileout, hwp_info, delimiter='       ', fmt=['%.0d','%.0d','%06.04f','%.0d'],header='FILE, HWP_POS, HWP_ANG, CYCLE')
    np.savetxt('index.txt', hwp_files, fmt='%d')
    print('HWP info saved to '+hwpfileout)

    #print('work in progress, no hwp log for you yet')

import sys
if __name__ == '__main__':
    if len(sys.argv)<2:
        datadir = input('Please specify file numbers for which you want to get HWP angles \n')
    else:
        datadir = sys.argv[1]
    hwp_match(datadir)