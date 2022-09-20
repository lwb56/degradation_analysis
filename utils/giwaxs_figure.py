import sys, os
import glob
import matplotlib.pylab as plt
import numpy as np
import scipy as sp
import pandas as pd
import PIL as pl
from imageio import imread


# from scipy.misc import imread


def plot_giwaxs(SciAnalysis_PATH, cal, mask, image_file, vmin, vmax):
    '''
    This function is for making the plotting for the GIWAXS image in python

    INPUTS:
        SciAnalysis_Path - the path where the SciAnalysis package is downloaded
        cal - a SciAnalysis object where the GIWAXS experimental parameters are kept
        mask - a SciAnalysis object of the mask file used in the experiment
        image_file - the GIWAXS pattern of interest
    Outputs:
        GIWAXS_fig - a matplotlib figure object which contain the GIWAXS image in q space
    '''
    SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)
    
    # importing the SciAnalysis package
    from SciAnalysis import tools
    from matplotlib import colors
    from matplotlib import ticker
    from matplotlib.ticker import AutoMinorLocator

    import matplotlib as mpl
    mpl.rcParams["axes.linewidth"] = 5
    mpl.rcParams["axes.labelsize"] = 40
    mpl.rcParams["axes.labelweight"] = 'bold'
    mpl.rcParams['xtick.major.size'] = 14
    mpl.rcParams['xtick.minor.size'] = 10
    mpl.rcParams['ytick.major.size'] = 14
    mpl.rcParams['ytick.minor.size'] = 10
    mpl.rcParams['xtick.major.width'] = 5
    mpl.rcParams['ytick.major.width'] = 5
    mpl.rcParams['xtick.minor.width'] = 4
    mpl.rcParams['ytick.minor.width'] = 4
    mpl.rcParams['font.size'] = 38
    mpl.rcParams["font.family"] = 'Arial'
    mpl.rcParams["font.weight"] = 'normal'
    mpl.rcParams['text.usetex'] = False


    # loading the image using PIL
    data = pl.Image.open(image_file)

    # image = imread(image_file)

    cal._generate_qxyz_maps()

    image = np.array(data)

    ##### Normalizing the image and then perform log
    # image = image/np.max(image)*100

    giwaxs = np.log(np.abs(image+11))
    giwaxs = giwaxs**2

    # giwaxs = image

    # print(image.data)
    # using the calibration file, we can extract the q maps specifically (Qr vs QZ)
    fig = plt.figure(figsize = (10,9))
    ax = fig.gca()
    # pc = ax.pcolor(cal.qr_map_data, cal.qz_map_data,giwaxs, cmap = 'plasma', norm= colors.Normalize(vmin= None, vmax = None))
    
    pc = ax.pcolor(cal.qr_map_data, cal.qz_map_data,giwaxs, cmap = 'jet', vmin =vmin, vmax = vmax)

    ax.set_xlim([-.15,2.2])
    ax.set_ylim([-.1,2.2])

    # ax.set_xlim([-2.2,2.2])
    # ax.set_ylim([-.1,2.2])


    ax.set_yticks(np.arange(0,2.2,.5))
    ax.set_xticks(np.arange(0,2.2,.5))
    # fig.colorbar(pc, ax = ax)
    # ax.tick_params(axis='x', which='minor', bottom=True)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel("$q_r (\AA^{-1})$")
    ax.set_ylabel("$q_z (\AA^{-1})$")
    fig.colorbar(pc)


    return None

