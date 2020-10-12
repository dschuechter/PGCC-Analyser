import random
import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import healpy as hp
from scipy import odr
from scipy.optimize import curve_fit
from scipy import stats
import pandas as pd

import sys
#Global Parameters
nside = 1024
fwhm=np.radians(0.27)
name = "none"
global_labels = ["NHI / cm$^{-2}$", "$I_{3THz}$ / MJy/sr", "$I_{857GHz}$ / MJycm$^2$/sr", "$I_{545GHz}$ / MJycm$^2$/sr","$I_{353GHz}$ /  MJycm$^2$/sr","PGCC", "Color $C$", "$T_{Dust, R1}$ / K", "$T_{Dust, R2}$ / K", "$W_{CO}$ / K $\cdot$ km/s"] 
global_labels_no_pgcc = ["NHI / cm$^{-2}$", "$I_{3THz}$ / MJycm$^2$/sr", "$I_{857GHz}$ / MJycm$^2$/sr", "$I_{545GHz}$ / MJycm$^2$/sr","$I_{353GHz}$ /  MJycm$^2$/sr", "Color $C$", "$T_{Dust, R1}$ / K", "$T_{Dust, R2}$ / K", "$W_{CO}$ / K $\cdot$ km/s"] 
global_names = ["Wasserstoffsäulendichte NHI", "Intensität $I_{3THz}$", "Intensität $I_{857GHz}$", "Intensität $I_{545GHz}$", "Intensität $I_{353GHz}$", "Color $C$", "Staubtemperatur $T_{Dust, R1}$", "Staubtemperatur $T_{Dust, R2}$", "Helligkeitstemperatur $W_{CO}$"] 
filetype = ".pdf"

# Debugging Tools
def gen_test_map(data,intensity,debugging, maximum):
    if(debugging==True):
        test_glon = data["Longitude"].values
        test_glat = -data["Latitude"].values+90
        test_intensity = data[intensity]
        pixel_indices = hp.ang2pix(nside, np.radians(test_glat), np.radians(test_glon))
        test_map = np.zeros(hp.nside2npix(nside))
        test_map[pixel_indices] = test_intensity
        if maximum!=0:
            hp.mollview(test_map, max=maximum, title="", cbar=False)
        else:
            hp.mollview(test_map, title="", cbar=False)
def debug_print(string,var, debugging):
    if debugging==True: print(string,":",var)
