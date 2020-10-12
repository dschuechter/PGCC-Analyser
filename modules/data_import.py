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
import os

#Global Parameters
nside = 1024
fwhm=np.radians(0.27)
name = "none"
global_labels = ["NHI / cm$^{-2}$", "$I_{3THz}$ / MJy/sr", "$I_{857GHz}$ / MJycm$^2$/sr", "$I_{545GHz}$ / MJycm$^2$/sr","$I_{353GHz}$ /  MJycm$^2$/sr","PGCC", "Color $C$", "$T_{Dust, R1}$ / K", "$T_{Dust, R2}$ / K", "$W_{CO}$ / K $\cdot$ km/s"]
global_labels_no_pgcc = ["NHI / cm$^{-2}$", "$I_{3THz}$ / MJycm$^2$/sr", "$I_{857GHz}$ / MJycm$^2$/sr", "$I_{545GHz}$ / MJycm$^2$/sr","$I_{353GHz}$ /  MJycm$^2$/sr", "Color $C$", "$T_{Dust, R1}$ / K", "$T_{Dust, R2}$ / K", "$W_{CO}$ / K $\cdot$ km/s"]
global_names = ["Wasserstoffsäulendichte NHI", "Intensität $I_{3THz}$", "Intensität $I_{857GHz}$", "Intensität $I_{545GHz}$", "Intensität $I_{353GHz}$", "Color $C$", "Staubtemperatur $T_{Dust, R1}$", "Staubtemperatur $T_{Dust, R2}$", "Helligkeitstemperatur $W_{CO}$"]
filetype = ".pdf"


def smooth(name, sigma):
    return hp.smoothing(name, sigma=(fwhm**2-sigma**2))
def downgrade(name, sigma):
    tmp = hp.ud_grade(name, nside_out=nside)
    return smooth(tmp, sigma)

def adjust_files():

    #PLANCK_353
    PLANCK_353 = hp.read_map("data/HFI_SkyMap_353_2048_R3.01_full.fits")
    PLANCK_353 = downgrade(PLANCK_353, 0.0012828170002158253)

    #PLANCK_545
    PLANCK_545 = hp.read_map("data/HFI_SkyMap_545_2048_R3.01_full.fits")
    PLANCK_545 = downgrade(PLANCK_545, 0.0013002702927357684)

    #PLANCK_857
    PLANCK_857 = hp.read_map("data/HFI_SkyMap_857_2048_R3.01_full.fits")
    PLANCK_857 = downgrade(PLANCK_857, 0.0012304571226559957)

    #IRAS_3000
    #NaN Values zerstören das smoothing der Daten. Deswegen werden diese noch konvertiert. Sie liegen nicht in der
    #beobachteten Draco Region
    IRIS_3000 = np.nan_to_num(hp.read_map("data/IRIS_nohole_4_1024_v2.fits"))
    IRIS_3000 = smooth(IRIS_3000, 0.001250819297262596)

    #Dust_Temp Revision 1 http://pla.esac.esa.int/pla-sl/data-action?MAP.MAP_OID=9300
    Dust_Temp_R1 = hp.read_map("data/HFI_CompMap_ThermalDustModel_2048_R1.20.fits", field=4)
    Dust_Temp_R1 = downgrade(Dust_Temp_R1, 0.0014544410433286001)

    #Dust_Temp http://pla.esac.esa.int/pla/#results searchterm "temperature"
    Dust_Temp_R2 = hp.read_map("data/COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.00.fits")
    Dust_Temp_R2 = downgrade(Dust_Temp_R2, 0.0014544410433286001)

    #Dame CO Map https://lambda.gsfc.nasa.gov/product/foreground/fg_dameco_map.cfm
    Dame = hp.read_map("data/lambda_wco_dht2001.fits")
    Dame = downgrade(Dame, np.radians(0.125))

    #Write adjusted Datasets
    fits.writeto("data/Comprimised/PLANCK_353.fits", PLANCK_353, overwrite=True)
    fits.writeto("data/Comprimised/PLANCK_545.fits", PLANCK_545, overwrite=True)
    fits.writeto("data/Comprimised/PLANCK_857.fits", PLANCK_857, overwrite=True)
    fits.writeto("data/Comprimised/IRIS_3000.fits", IRIS_3000, overwrite=True)
    fits.writeto("data/Comprimised/Dust_Temp_R1.fits", Dust_Temp_R1, overwrite=True)
    fits.writeto("data/Comprimised/Dust_Temp_R2.fits", Dust_Temp_R2, overwrite=True)
    fits.writeto("data/Comprimised/Dame.fits", Dame, overwrite=True)


    return PLANCK_353, PLANCK_545, PLANCK_857, IRIS_3000, Dust_Temp_R1, Dust_Temp_R2, Dame



def import_files():
    #PLANCK_353
    PLANCK_353 = fits.getdata("data/Comprimised/PLANCK_353.fits")

    #PLANCK_545
    PLANCK_545 = fits.getdata("data/Comprimised/PLANCK_545.fits")

    #PLANCK_857
    PLANCK_857 = fits.getdata("data/Comprimised/PLANCK_857.fits")

    #IRAS_3000
    IRIS_3000 = fits.getdata("data/Comprimised/IRIS_3000.fits")

    #Dust_Temp_R1
    Dust_Temp_R1 = fits.getdata("data/Comprimised/Dust_Temp_R1.fits")

    #Dust_Temp_R2
    Dust_Temp_R2 = fits.getdata("data/Comprimised/Dust_Temp_R2.fits")


    #CO
    Dame = fits.getdata("data/Comprimised/Dame.fits")

    return PLANCK_353, PLANCK_545, PLANCK_857, IRIS_3000, Dust_Temp_R1, Dust_Temp_R2, Dame

def import_data():
    if os.path.exists(os.getcwd()+"/data/comprimised"):
        PLANCK_353, PLANCK_545, PLANCK_857, IRIS_3000, Dust_Temp_R1, Dust_Temp_R2, Dame = import_files()
    else:
        os.mkdir(os.getcwd()+"/data/comprimised")
        PLANCK_353, PLANCK_545, PLANCK_857, IRIS_3000, Dust_Temp_R1, Dust_Temp_R2, Dame = adjust_files()
    #PGCC
    PGCC_in = fits.getdata('data/browse_results.fits', 1)
    PGCC_data = Table(PGCC_in)
    PGCC_glon = PGCC_data["LII"]
    PGCC_glat = -PGCC_data["BII"]
    PGCC_glat += 90
    PGCC_Intensity = np.ones(len(PGCC_glat))
    pixel_indices = hp.ang2pix(nside, np.radians(PGCC_glat), np.radians(PGCC_glon))
    PGCC = np.zeros(hp.nside2npix(nside))
    PGCC[pixel_indices] = PGCC_Intensity

    #HI4PI
    HI4PI_in = fits.getdata('data/NHI_HPX.fits', 1)
    HI4PI_data = Table(HI4PI_in)
    HI4PI_glon = HI4PI_data["GLON"]
    HI4PI_glat = -HI4PI_data["GLAT"]
    HI4PI_glat += 90
    HI4PI_Intensity = HI4PI_in["NHI"]
    pixel_indices = hp.ang2pix(nside, np.radians(HI4PI_glat), np.radians(HI4PI_glon))
    HI4PI = np.zeros(hp.nside2npix(nside))
    HI4PI[pixel_indices] = HI4PI_Intensity
    #Generiere Globale Koordinaten, s.d. diese nach der verwendung in Cartview ebenfalls erhalten bleiben
    NPIX = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(NPIX))
    global_lon = np.rad2deg(phi).astype(np.float64)
    global_lat = (90. - np.rad2deg(theta)).astype(np.float64)

    Farbkarte = (PLANCK_857-IRIS_3000)/(IRIS_3000+PLANCK_857)
    all_sky_table = {"Longitude":global_lon, "Latitude":global_lat,"NHI": HI4PI.ravel(), "IRIS_3000GHz":IRIS_3000.ravel(), "PLANCK_857GHz":PLANCK_857.ravel(), "PLANCK_545GHz":PLANCK_545.ravel(), "PLANCK_353GHz":PLANCK_353.ravel(),"PGCC":PGCC.ravel(), "Color":Farbkarte.ravel(), "Dust_Temp_R1":Dust_Temp_R1, "Dust_Temp_R2":Dust_Temp_R2, "CO":Dame}
    all_sky_table = pd.DataFrame(data=all_sky_table)
    all_sky_PGCC_Table = all_sky_table.loc[all_sky_table['PGCC'] == 1]
    all_sky_PGCC_Table = all_sky_PGCC_Table.drop_duplicates()

    return HI4PI,PLANCK_353, PLANCK_545, PLANCK_857, IRIS_3000, Dust_Temp_R1, Dust_Temp_R2, Dame, all_sky_table
