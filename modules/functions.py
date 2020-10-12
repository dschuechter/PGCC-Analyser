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
import modules.debugging as debugging_tools
import os
import shutil
#Global Parameters
nside = 1024
fwhm=np.radians(0.27)
name = "none"
global_labels = ["NHI / cm$^{-2}$", "$I_{3THz}$ / MJy/sr", "$I_{857GHz}$ / MJycm$^2$/sr", "$I_{545GHz}$ / MJycm$^2$/sr","$I_{353GHz}$ /  MJycm$^2$/sr","PGCC", "Color $C$", "$T_{Dust, R1}$ / K", "$T_{Dust, R2}$ / K", "$W_{CO}$ / K $\cdot$ km/s"]
global_labels_no_pgcc = ["NHI / cm$^{-2}$", "$I_{3THz}$ / MJycm$^2$/sr", "$I_{857GHz}$ / MJycm$^2$/sr", "$I_{545GHz}$ / MJycm$^2$/sr","$I_{353GHz}$ /  MJycm$^2$/sr", "Color $C$", "$T_{Dust, R1}$ / K", "$T_{Dust, R2}$ / K", "$W_{CO}$ / K $\cdot$ km/s"]
global_units = ["cm$^{-2}$", "MJy/sr", "MJy/sr", "MJy/sr", "K$_{CMB}$", "", "K", "K", "K $\cdot$ km/s"]
global_names = ["Hydrogencolumndensity NHI", "Intensity $I_{3THz}$", "Intensity $I_{857GHz}$", "Intensity $I_{545GHz}$", "Intensity $I_{353GHz}$", "Color $C$", "Dusttemperature $T_{Dust, R1}$", "Dusttemperature $T_{Dust, R2}$", "Brightnesstemperature $W_{CO}$"]
filetype = ".pdf"
filepath = "./output/"


#Funktion um Daten aus gewünschten Himmelsregionen zu extrahieren
def data_selector(all_sky_table,lon_lim_min, lon_lim_max, lat_lim_min, lat_lim_max, inverted, lon_steps, lat_steps, debugging):
    if(inverted == True):
        selected_data = all_sky_table[(all_sky_table["Longitude"]>=lon_lim_min) & (all_sky_table["Longitude"]<=lon_lim_max) & ((all_sky_table["Latitude"]<=lat_lim_min) | (all_sky_table["Latitude"]>=lat_lim_max))]
    else:
        selected_data = all_sky_table[(all_sky_table["Longitude"]>=lon_lim_min) & (all_sky_table["Longitude"]<=lon_lim_max) & (all_sky_table["Latitude"]>=lat_lim_min)&(all_sky_table["Latitude"]<=lat_lim_max)]
    data_lat_slices = []
    if(lat_steps != 0):
        if inverted == True:
            loop_time = (180-np.abs(lat_lim_min)-np.abs(lat_lim_max))/lat_steps
            current_loop = 0
            debugging_tools.debug_print("loop_time_lat",loop_time, debugging)
            for i in range(0, int(loop_time)):
                data_lat_slices.append([])
                if(lat_lim_max+(i+1)*lat_steps<=90):
                    data_lat_slices[i] = selected_data[(selected_data["Latitude"]>=lat_lim_max+i*lat_steps) &
                                                  (selected_data["Latitude"]<=lat_lim_max+(i+1)*lat_steps) &
                                                   (lat_lim_max+(i+1)*lat_steps<=90)]
                    current_loop += 1
                else:
                    j = i-current_loop
                    data_lat_slices[i] = selected_data[(selected_data["Latitude"]<=lat_lim_min-j*lat_steps)&
                                                    (selected_data["Latitude"]>=lat_lim_min-(j+1)*lat_steps) &
                                                   (lat_lim_min-(j+1)*lat_steps>=-90)]

                #gen_test_map(data_lat_slices[i], "NHI", debugging)

        else:
            loop_time = (lat_lim_max-lat_lim_min)/lat_steps
            current_loop = 0
            debugging_tools.debug_print("loop_time_lat",loop_time, debugging)
            for i in range(0, int(loop_time)):
                data_lat_slices.append([])
                data_lat_slices[i] = selected_data[(selected_data["Latitude"]>=lat_lim_min+i*lat_steps) &
                                              (selected_data["Latitude"]<=lat_lim_min+(i+1)*lat_steps)]
                #gen_test_map(data_lat_slices[i], "NHI", debugging)
    data_lon_slices = []
    if(lon_steps != 0):
        if(lat_steps != 0):
            loop_time = (lon_lim_max-lon_lim_min)/lon_steps
            debugging_tools.debug_print("lon_slicing_loop_time", loop_time, debugging)
            current_loop = 0
            for i in range(0, len(data_lat_slices)):
                for j in range(0, int(loop_time)):
                    data_lon_slices.append([])
                    data_lon_slices[current_loop] = data_lat_slices[i][(data_lat_slices[i]["Longitude"]>=lon_lim_min+j*lon_steps)&
                                                             (data_lat_slices[i]["Longitude"]<=lon_lim_min+(j+1)*lon_steps)
                                                             ]
                    current_loop +=1
            #for z in range(0, len(data_lon_slices)): gen_test_map(data_lon_slices[z], "NHI", debugging)

        else:
            loop_time = (lon_lim_max-lon_lim_min)/lon_steps
            current_loop = 0
            debugging_tools.debug_print("loop_time_lat",loop_time, debugging)
            for i in range(0, int(loop_time)):
                data_lon_slices.append([])
                data_lon_slices[i] = selected_data[(selected_data["Longitude"]>=lon_lim_min+i*lon_steps) &
                                                  (selected_data["Longitude"]<=lon_lim_min+(i+1)*lon_steps)]
                #gen_test_map(data_lon_slices[i], "NHI", debugging)
    if(lon_steps == 0 and lat_steps == 0):
        data_slices = selected_data
    if(len(data_lon_slices)!=0):
        data_slices = selected_data
    if(len(data_lat_slices)!=0 and len(data_lon_slices)==0):
        data_slices = data_lat_slices
    return selected_data, data_slices

#Funktion um aus bestehenden Datensatz Daten zu filtern
def filter_data(data, filter_option):
    filter_lon = [0,0]
    filter_lat = [0,0]
    if filter_option == "LMC":
        filter_lon[0]=273.94
        filter_lon[1]=283.12
        filter_lat[0]=-37.69
        filter_lat[1]=-29.51
    if filter_option == "SMC":
        filter_lon[0]=300.29
        filter_lon[1]=305.05
        filter_lat[0]=-46.86
        filter_lat[1]=-42.08
    if filter_option == "orion_nebula":
        filter_lon[0]=205
        filter_lon[1]=215
        filter_lat[0]=-27
        filter_lat[1]=-15
    if filter_option == "orion_super_bubble":
        filter_lon[0]=207.5
        filter_lon[1]=210
        filter_lat[0]=-21
        filter_lat[1]=-18
    if filter_option == "all":
        filter_lon[0]=207.5
        filter_lon[1]=210
        filter_lat[0]=-21
        filter_lat[1]=-18
        negative_map, negative_map_sliced = data_selector(data,filter_lon[0], filter_lon[1], filter_lat[0], filter_lat[1], False, 0, 0, False)
        data = pd.concat([data, negative_map, negative_map]).drop_duplicates(keep=False)
        filter_lon[0]=300.29
        filter_lon[1]=305.05
        filter_lat[0]=-46.86
        filter_lat[1]=-42.08
        negative_map, negative_map_sliced = data_selector(data,filter_lon[0], filter_lon[1], filter_lat[0], filter_lat[1], False, 0, 0, False)
        data = pd.concat([data, negative_map, negative_map]).drop_duplicates(keep=False)
        filter_lon[0]=273.94
        filter_lon[1]=283.12
        filter_lat[0]=-37.69
        filter_lat[1]=-29.51
        negative_map, negative_map_sliced = data_selector(data,filter_lon[0], filter_lon[1], filter_lat[0], filter_lat[1], False, 0, 0, False)
        data = pd.concat([data, negative_map, negative_map]).drop_duplicates(keep=False)
    if filter_option != "MC":
        negative_map, negative_map_sliced = data_selector(data,filter_lon[0], filter_lon[1], filter_lat[0], filter_lat[1], False, 0, 0, False)
        data = pd.concat([data, negative_map, negative_map]).drop_duplicates(keep=False)
        test_map = debugging_tools.gen_test_map(data, "NHI",True, 0)
    else:
        filter_lon[0]=273.94
        filter_lon[1]=283.12
        filter_lat[0]=-37.69
        filter_lat[1]=-29.51
        negative_map, negative_map_sliced = data_selector(data,filter_lon[0], filter_lon[1], filter_lat[0], filter_lat[1], False, 0, 0, False)
        data = pd.concat([data, negative_map, negative_map]).drop_duplicates(keep=False)
        filter_lon[0]=300.29
        filter_lon[1]=305.05
        filter_lat[0]=-46.86
        filter_lat[1]=-42.08
        negative_map, negative_map_sliced = data_selector(data,filter_lon[0], filter_lon[1], filter_lat[0], filter_lat[1], False, 0, 0, False)
        data = pd.concat([data, negative_map, negative_map]).drop_duplicates(keep=False)
        test_map = debugging_tools.gen_test_map(data, "NHI",True, 0)
    return data




#Funktion um einen simplen scatter plot zu erstellen
def scatter_plot(x, y,x_label, y_label, save_name, lat_steps,lon_steps, show_figure):
    plt.figure()
    if(lon_steps!=0 or lat_steps!=0):
        for i in range(0, len(x)):
            plt.scatter(x[i],y[i], color ="black")

    else:
        plt.scatter(x,y, color = "black")
    if(show_figure==True):
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.savefig(filepath+name+"_"+save_name+"_"+"scatter"+filetype,bbox_inches='tight')
        plt.show()

#Funtion um einen linearen Fit zu erstelleng
def lin_fit(xdata,ydata):
    z = np.polyfit(xdata, ydata, 1)
    f = np.poly1d(z)
    corrcoeff = np.corrcoef(xdata, ydata)[0,1]
    return f[1], f[0], corrcoeff
#Funktion um einen linearen Fit zu erstellen und Fehler anzugeben
def func_lin(p,x):
    a,b=p
    return a*x+b
def anpassung_xerr(function, x, y, x_error, presets, plot,customlabel):
    model = odr.Model(function)
    data = odr.RealData(x, y, sx=x_error)
    out = odr.ODR(data, model, beta0=presets).run()

    popt = out.beta
    perr = out.sd_beta

    if plot == True:
        x_fit = np.linspace(min(x), max(x), 10000)
        y_fit = function(popt, x_fit)
        if customlabel == False:
            plt.plot(x_fit, y_fit)
        else:
            plt.plot(x_fit, y_fit,label=customlabel)
    return popt,perr

def anpassung_yerr(function, x, y, y_error, presets, plot,customlabel):
    model = odr.Model(function)
    data = odr.RealData(x, y, sy=y_error)
    out = odr.ODR(data, model, beta0=presets).run()

    popt = out.beta
    perr = out.sd_beta

    if plot == True:
        x_fit = np.linspace(min(x), max(x), 10000)
        y_fit = function(popt, x_fit)
        if customlabel == False:
            plt.plot(x_fit, y_fit)
        else:
            plt.plot(x_fit, y_fit,label=customlabel)
    return popt,perr

def anpassung_no_err(function, x, y, presets, plot,customlabel):
    model = odr.Model(function)
    data = odr.RealData(x, y)
    out = odr.ODR(data, model, beta0=presets).run()

    popt = out.beta
    perr = out.sd_beta

    return popt,perr

#Funktion um übergebene Daten zu korrelieren
def korrelation(data, x_header, y_header, fit_type, x_label, y_label, show_figure):
    m = 0
    b = 0
    c = 0
    if(fit_type == "lin_fit"):
        m,b,c = lin_fit(data[x_header], data[y_header])
        if(show_figure==True):
            plt.figure()
            x_fit = np.linspace(np.amin(data[x_header]), np.amax(data[x_header]), 1000)
            plt.scatter(data[x_header], data[y_header], s=5)
            y_fit = m*x_fit+b
            plt.plot(x_fit, y_fit, "-", color="black", label="$f(x)=%s\cdot x + %s$\nKorrelationskoeffizient: %.3f\nKoord: $Lon:[%.3f,%.3f], Lat: [%.3f,%.3f]$"%(format(m, '.2e'),format(b, '.2e'),c, np.amin(data["Longitude"]),np.amax(data["Longitude"]),np.amin(data["Latitude"]),np.amax(data["Latitude"])))
            for i in range(0, len(global_labels)):
                if x_label in global_labels[i]:
                    x_label = global_labels[i]
                if y_label in global_labels[i]:
                    y_label = global_labels[i]
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.grid()
            plt.legend(loc="upper right")
            plt.savefig(filepath+name+"_"+x_header+"_"+y_header+"_"+"korrelation.png")
            plt.show()
    return m,b,np.abs(c)

def TWO_D_korrelation(data, x_header, y_header, z_header, fit_type, x_label, y_label, z_label):
    m = 0
    b = 0
    c = 0
    if(fit_type == "lin_fit"):
        m,b,c = lin_fit(data[x_header], data[y_header])
        plt.figure()
        x_fit = np.linspace(np.amin(data[x_header]), np.amax(data[x_header]), 1000)
        plt.scatter(data[x_header], data[y_header],c=data[z_header], s=5)
        y_fit = m*x_fit+b
        plt.plot(x_fit, y_fit, "-", color="black", label="$f(x)=%.3f\cdot x + %.3f$\nKorrelationskoeffizient: %.3f\nKoord: $Lon:[%.3f,%.3f], Lat: [%.3f,%.3f]$"%(m,b,c, np.amin(data["Longitude"]),np.amax(data["Longitude"]),np.amin(data["Latitude"]),np.amax(data["Latitude"])))
        for i in range(0, len(global_labels)):
            if x_label in global_labels[i]:
                x_label = global_labels[i]
            if y_label in global_labels[i]:
                y_label = global_labels[i]
            if z_label in global_labels[i]:
                z_label = global_labels[i]
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.legend(loc="upper right")
        plt.grid()
        plt.colorbar(label=z_label)
        plt.savefig(filepath+name+"_"+x_header+"_"+y_header+"_"+"korrelation.png")
        plt.show()
    return m,b,np.abs(c)
#Funktion um eingegebene Daten zu mitteln
def mitteln(data, header_name, lon_steps, lat_steps):
    avrg = []
    if(lon_steps!=0 or lat_steps!=0):
        for i in range(0,len(data)):
            avrg.append([])
            avrg[i] = np.average(data[i][header_name])
    else:
        avrg = np.average(data[header_name])
    return avrg

#Funktion um einen ausgewählten Himmelsbereich für die verschiedenen Karteninhalte darzustellen
def show_pics(data,lon_steps, lat_steps, pairs, no_pgcc_list):
    if(lon_steps==0 and lat_steps==0):
        fig, axes = plt.subplots(nrows=3, ncols=3)
        fig.set_figheight(10)
        fig.set_figwidth(15)
        z = 2
        head = data.columns.values
        for i in range(0,3):
            for j in range(0,3):
                if "PGCC" in head[z]:
                    z += 1
                if "Dust_Temp" in head[z]:
                    img = axes[i,j].scatter(x=data["Longitude"],y=data["Latitude"],c=data[head[z]],vmin=10, vmax=25,s=1, cmap='viridis')
                    #data.plot.scatter(ax=axes[i,j],x="Longitude",y="Latitude",c=head[z],vmin=10, vmax=25,s=1,subplots=True, colormap='viridis',antialiaseds=True)
                elif("NHI" in head[z]):
                    img = axes[i,j].scatter(x=data["Longitude"],y=data["Latitude"],c=data[head[z]],s=1, cmap='viridis')
                    img = axes[i,j].scatter(x=data[data["PGCC"]==1]["Longitude"],y=data[data["PGCC"]==1]["Latitude"],c="white",s=4)
                    if pairs != False:
                        pairs = [ item for elem in pairs for item in elem]
                        img = axes[i,j].scatter(x=data[data["PGCC"]==1].iloc[pairs]["Longitude"],y=data[data["PGCC"]==1].iloc[pairs]["Latitude"],c="black",s=4)
                    if no_pgcc_list != False:
                        #no_pgcc_list = [ item for elem in no_pgcc_list for item in elem]
                        img = axes[i,j].scatter(x=data[data["PGCC"]==1].iloc[no_pgcc_list]["Longitude"],y=data[data["PGCC"]==1].iloc[no_pgcc_list]["Latitude"],c="red",s=4)
                    #data.plot.scatter(ax=axes[i,j],x="Longitude",y="Latitude",c=head[z],s=1,subplots=True, colormap='viridis')
                    #data[data["PGCC"]==1].plot.scatter(ax=axes[i,j],x="Longitude",y="Latitude",c="red",s=1)
                else:
                    #data.plot.scatter(ax=axes[i,j],x="Longitude",y="Latitude",c=head[z],s=1,subplots=True, colormap='viridis')
                    img = axes[i,j].scatter(x=data["Longitude"],y=data["Latitude"],c=data[head[z]],s=1, cmap='viridis')
                axes[i,j].set_xlabel("Longitude")
                axes[i,j].set_ylabel("Latitude")
                plt.colorbar(img, ax=axes[i,j], label=global_labels[z-2])
                z += 1
        fig.tight_layout()
        plt.savefig(filepath+name+"_"+"pics.png")
        plt.show()

def Remove_Duplicates(lst):
     return ([list(i) for i in {*[tuple(sorted(i)) for i in lst]}])
###
# Werkzeuge um die PGCC's zu untersuchen
###
def PGCC_eight_pix(all_sky_Table,lon, lat):
        color_avrg = []
        NHI_avrg = []
        IRIS_3000_avrg = []
        PLANCK_353_avrg = []
        PLANCK_545_avrg = []
        PLANCK_857_avrg = []
        dust_temp_R_1_avrg = []
        dust_temp_R_2_avrg = []
        CO_avrg = []
        for i in range(0,len(lon)):
            coords = hp.get_all_neighbours(nside, lon.values[i], lat.values[i], lonlat=True)
            NHI_avrg.append(np.average(all_sky_Table.iloc[coords]["NHI"]))
            color_avrg.append(np.average(all_sky_Table.loc[coords]["Color"]))
            dust_temp_R_1_avrg.append(np.average(all_sky_Table.loc[coords]["Dust_Temp_R1"]))
            dust_temp_R_2_avrg.append(np.average(all_sky_Table.loc[coords]["Dust_Temp_R2"]))
            IRIS_3000_avrg.append(np.average(all_sky_Table.loc[coords]["IRIS_3000GHz"]))
            PLANCK_353_avrg.append(np.average(all_sky_Table.loc[coords]["PLANCK_353GHz"]))
            PLANCK_545_avrg.append(np.average(all_sky_Table.loc[coords]["PLANCK_545GHz"]))
            PLANCK_857_avrg.append(np.average(all_sky_Table.loc[coords]["PLANCK_857GHz"]))
            CO_avrg.append(np.average(all_sky_Table.loc[coords]["CO"]))

        table = {"Longitude":lon, "Latitude":lat, "NHI_avrg": NHI_avrg, "IRIS_3000_avrg":IRIS_3000_avrg,"PLANCK_857_avrg":PLANCK_857_avrg,"PLANCK_545_avrg":PLANCK_545_avrg,"PLANCK_353_avrg":PLANCK_353_avrg,"Color_avrg": color_avrg, "Dust_Temp_R_1_avrg":dust_temp_R_1_avrg, "Dust_Temp_R_2_avrg":dust_temp_R_2_avrg, "CO_avrg": CO_avrg}
        return pd.DataFrame(data=table)

def wideview(lon, lat, c, s, label):
    for i in range(0,len(lon)):
        if(lon.iloc[i]>180):
            lon.iloc[i]=lon.iloc[i]-360
    if label != False:
        plt.scatter(lon, lat, c=c, s=s, label=label)
    else:
        plt.scatter(lon, lat,c=c, s=s)

def find_PGCCs_nearby(data, dist,Color_tollerance):
    pairs = []
    pairs_sorted = []
    skip = []
    sorted_index=0
    for i in range(0, len(data)):
        pairs.append([])
        #Vorzeichen Problematik!!!!!!
        tmp_lon = data["Longitude"].values-data["Longitude"].values[i]
        tmp_lat = data["Latitude"].values-data["Latitude"].values[i]
        for j in range(0, len(data)):
            if(np.sqrt(np.abs(tmp_lon[j])**2+np.abs(tmp_lat[j])**2)<=dist and np.abs(tmp_lon[j])+np.abs(tmp_lat[j])!=0):
                pairs[i].append(j)
        pairs[i].append(i)

    ##
    # Entferne nicht vorhandene Paare
    ##
    pairs = [x for x in pairs if len(x) > 1]

    ##
    # Entferne Duplikate
    ##
    pairs_sorted = Remove_Duplicates(pairs)

    print_pairs = [item for elem in pairs_sorted for item in elem]
    print_pairs = list(dict.fromkeys(print_pairs))
    print("There are ", len(print_pairs) , " PGCC's in a neighbourhood")

    if len(pairs)!= 0:
        ###
        # Prüfe ob die Paare die gleiche Farbe aufweisen -> Wahrscheinlich ein und das selbe Objekt
        ###
        Color_pairs = []
        for i in range(0, len(pairs_sorted)):
            tmp_pairs=[]
            for j in range(0, len(pairs_sorted[i])):
                for z in range(0, len(pairs_sorted[i])):
                    if pairs_sorted[i][z]!=pairs_sorted[i][j]:
                        if(data.iloc[pairs_sorted[i][j]]["Color_avrg"]<=Color_tollerance+data.iloc[pairs_sorted[i][z]]["Color_avrg"] and data.iloc[pairs_sorted[i][j]]["Color_avrg"]>=data.iloc[pairs_sorted[i][z]]["Color_avrg"]-Color_tollerance):
                            tmp_pairs.append(pairs_sorted[i][j])
                            tmp_pairs.append(pairs_sorted[i][z])
            tmp_pairs = list(dict.fromkeys(tmp_pairs))
            Color_pairs.append(tmp_pairs)
        Color_pairs = [x for x in Color_pairs if len(x) > 1]
        pairs_sorted_color = Color_pairs.copy()
        print_pairs = [item for elem in Color_pairs for item in elem]
        print_pairs = list(dict.fromkeys(print_pairs))
        print("Of all found PGCC's were ", len(print_pairs) , " matches in Color")


        ###
        # Find out how many objects are real PGCC's and/or part of big structures
        ###

        neighbour_list = print_pairs
        PGCC_objects = []
        objects_found = []
        count = 0
        ignore = []
        PGCC_objects = []
        object_count = 0

        for i in range(0, len(Color_pairs)):
            if i not in ignore:
                ignore.append(i)
                PGCC_objects.append([])
                PGCC_objects[object_count]=Color_pairs[i]
                for z in range(0,2):
                    for j in range(0, len(Color_pairs)):
                        if j not in ignore:
                            if((np.setdiff1d(Color_pairs[j],PGCC_objects[object_count])).tolist()!=Color_pairs[j]):
                                ignore.append(j)
                                PGCC_objects[object_count]=PGCC_objects[object_count]+(np.setdiff1d(Color_pairs[j],PGCC_objects[object_count])).tolist()
                                PGCC_objects[object_count]=list(dict.fromkeys(PGCC_objects[object_count]))
                object_count+=1

        print("")
        print("%d coherent structures have been found"%(len(PGCC_objects)))
        if len(PGCC_objects)< 20:
            print(PGCC_objects)
        durchschnittsquellenanzahl = 0
        maximalequellenanzahl = 0
        minimalequellenanzahl = len(PGCC_objects[0])
        neighbour_list_lengths = []
        for i in range(0, len(PGCC_objects)):
            durchschnittsquellenanzahl+=len(PGCC_objects[i])
            if len(PGCC_objects[i])<minimalequellenanzahl:
                minimalequellenanzahl = len(PGCC_objects[i])
            if len(PGCC_objects[i])>maximalequellenanzahl:
                maximalequellenanzahl = len(PGCC_objects[i])
            neighbour_list_lengths.append(len(PGCC_objects[i]))

        durchschnittsquellenanzahl=durchschnittsquellenanzahl/len(PGCC_objects)
        print("\nParameteres of the structures:\nAverage number of sources: %.2f\nMinimal amount of sources: %d\nMaxmimal amount of sources: %d\n"%(durchschnittsquellenanzahl, minimalequellenanzahl, maximalequellenanzahl))
        n_lenghts = [0]*(np.amax(neighbour_list_lengths)+1)
        for i in range(0, len(neighbour_list_lengths)):
            n_lenghts[neighbour_list_lengths[i]]+=1
        for i in range(0, len(n_lenghts)):
            if(n_lenghts[i]>0):
                print(n_lenghts[i]," structures with ", i, " sources found")
        plt.hist(x=neighbour_list_lengths, bins=100)
        plt.xlabel("Number of PGCC's in structures")
        plt.ylabel("Frequency")
        plt.savefig(filepath+name+"pgcc_neighbourhood_"+str(dist)+".png")
        plt.plot()

        expanded_pgcc_found = 0
        expanded_list = []
        for i in range(0, len(PGCC_objects)):
            tmp_lon = 0
            tmp_lat = 0
            for j in range(0, len(PGCC_objects[i])):
                tmp_lon += data.iloc[PGCC_objects[i][j]]["Longitude"]
                tmp_lat += data.iloc[PGCC_objects[i][j]]["Latitude"]
            expanded_pgcc = len(PGCC_objects[i])
            for j in range(0, len(PGCC_objects[i])):
                if(np.sqrt((data.iloc[PGCC_objects[i][j]]["Longitude"]-tmp_lon/len(PGCC_objects[i]))**2+(data.iloc[PGCC_objects[i][j]]["Latitude"]-tmp_lat/len(PGCC_objects[i]))**2<1.5)):
                    expanded_pgcc -= 1
            if  expanded_pgcc == 0:
                expanded_pgcc_found += 1
                expanded_list.append(PGCC_objects[i])
        expanded_pgcc = [item for elem in expanded_list for item in elem]
        expanded_pgcc = list(dict.fromkeys(expanded_pgcc))
        PGCC_objects_list = [item for elem in PGCC_objects for item in elem]
        no_pgcc_list = (np.setdiff1d(PGCC_objects_list,expanded_pgcc)).tolist()
        print("\nOf %d found PGCC's are likely %d real individual sources"%(len(data), len(data)-len(neighbour_list)))
        print("%d structures were found, that build up single PGCC's"%(expanded_pgcc_found))
        print("Probably %d listed PGCC's are no sources, but part of a bigger extended structures"%(len(no_pgcc_list)))
        if expanded_pgcc_found!=0:
            print(expanded_list)
            plt.figure(figsize=(21, 12))
            wideview(data["Longitude"], data["Latitude"], "blue", 0.1, "PGCC-Singlesources")
            wideview(data.iloc[expanded_list[0]]["Longitude"], data.iloc[expanded_list[0]]["Latitude"], "black", 0.3, "PGCC-Neighbourhoods")
            for i in range(1, len(expanded_list)):
                wideview(data.iloc[expanded_list[i]]["Longitude"], data.iloc[expanded_list[i]]["Latitude"], "black", 0.3, False)
            wideview(data.iloc[no_pgcc_list]["Longitude"], data.iloc[no_pgcc_list]["Latitude"], "red", 0.3, "Big structures")
            plt.legend()
            plt.xlabel("Longitude / Degree")
            plt.ylabel("Latitude / Degree")
            plt.savefig(filepath+name+"neighbourhood"+str(dist)+filetype,bbox_inches='tight')
            plt.show()
        return PGCC_objects,expanded_list
    else:
        print("No structures could be found!")
        return False,False


#####
# Method to find deviating PGCC's
#####
def pgcc_analyser(pgcc, pgcc_ref, show_figure):

    fig, axes = plt.subplots(nrows=3, ncols=3)
    fig.set_figheight(10)
    fig.set_figwidth(15)
    #axes[-1,-1].axis('off')
    mean = [0]*9
    std = [0]*9
    z = 2
    head = pgcc.columns.values
    for i in range(0,3):
        for j in range(0,3):
            mean[z-2], std[z-2] = pgcc_ref[head[z]].agg(['mean', 'std']).round(decimals=10)
            if global_units[z-2]!="":
                pgcc_ref[head[z]].plot.hist(ax=axes[i,j],grid=True, bins=100, label ="$\mu$=%s/%s\n$\sigma$=%s/%s"%(format(mean[z-2], '.2e'),global_units[z-2], format(std[z-2], '.2e'),global_units[z-2]),subplots=True)
            else:
                pgcc_ref[head[z]].plot.hist(ax=axes[i,j],grid=True, bins=100, label ="$\mu$=%s\n$\sigma$=%s"%(format(mean[z-2], '.2e'), format(std[z-2], '.2e')),subplots=True)
            axes[i,j].legend()
            axes[i,j].set_xlabel(global_labels_no_pgcc[z-2])
            axes[i,j].set_ylabel("Häufigkeit")
            axes[i,j].set_title(global_names[z-2])
            if show_figure == True:
                plt.grid(axis="x", alpha=0.75)
            z += 1
            if z == 11: break;
    fig.tight_layout()
    if show_figure == True:
        plt.savefig(filepath+name+"_"+"pgcc.png")
        plt.show()

    print("\n---Strange PGCC's---\n")
    #Find strange behaving PGCC's
    strange_PGCC = pgcc.copy()
    head = strange_PGCC.columns.values
    NHI_flag = [0]*len(strange_PGCC)
    IRIS_3000_flag = [0]*len(strange_PGCC)
    PLANCK_857_flag = [0]*len(strange_PGCC)
    PLANCK_545_flag = [0]*len(strange_PGCC)
    PLANCK_353_flag = [0]*len(strange_PGCC)
    Color_flag = [0]*len(strange_PGCC)
    Dust_R_1_flag = [0]*len(strange_PGCC)
    Dust_R_2_flag = [0]*len(strange_PGCC)
    CO_flag = [0]*len(strange_PGCC)
    sigma = 1
    for i in range(0,len(strange_PGCC)):
        if((strange_PGCC["NHI_avrg"].iloc[i]<=sigma*(mean[0]-std[0])) or (strange_PGCC["NHI_avrg"].iloc[i]>=sigma*(mean[0]+std[0]))):
            NHI_flag[i]=1
        if((strange_PGCC["IRIS_3000_avrg"].iloc[i]<=sigma*(mean[1]-std[1])) or (strange_PGCC["IRIS_3000_avrg"].iloc[i]>=sigma*(mean[1]+std[1]))):
            IRIS_3000_flag[i]=1
        if((strange_PGCC["PLANCK_857_avrg"].iloc[i]<=sigma*(mean[2]-std[2])) or (strange_PGCC["PLANCK_857_avrg"].iloc[i]>=sigma*(mean[2]+std[2]))):
            PLANCK_857_flag[i]=1
        if((strange_PGCC["PLANCK_545_avrg"].iloc[i]<=sigma*(mean[3]-std[3])) or (strange_PGCC["PLANCK_545_avrg"].iloc[i]>=sigma*(mean[3]+std[3]))):
             PLANCK_545_flag[i]=1
        if((strange_PGCC["PLANCK_353_avrg"].iloc[i]<=sigma*(mean[4]-std[4])) or (strange_PGCC["PLANCK_353_avrg"].iloc[i]>=sigma*(mean[4]+std[4]))):
             PLANCK_353_flag[i]=1
        if((strange_PGCC["Color_avrg"].iloc[i]<=sigma*(mean[5]-std[5])) or (strange_PGCC["Color_avrg"].iloc[i]>=sigma*(mean[5]+std[5]))):
            Color_flag[i]=1
        if(strange_PGCC["Dust_Temp_R_1_avrg"].iloc[i]>=sigma*(mean[6]+std[6])):
            Dust_R_1_flag[i]=1
        if(strange_PGCC["Dust_Temp_R_2_avrg"].iloc[i]>=sigma*(mean[7]+std[7])):
            Dust_R_2_flag[i]=1
        if(strange_PGCC["CO_avrg"].iloc[i]>=sigma*(mean[8]+std[8])):
            CO_flag[i]=1
    strange_PGCC["NHI_flag"]=NHI_flag
    strange_PGCC["IRIS_3000_flag"]=IRIS_3000_flag
    strange_PGCC["PLANCK_857_flag"]=PLANCK_857_flag
    strange_PGCC["PLANCK_545_flag"]=PLANCK_545_flag
    strange_PGCC["PLANCK_353_flag"]=PLANCK_353_flag
    strange_PGCC["Color_flag"]=Color_flag
    strange_PGCC["Dust_Temp_R_1_flag"]=Dust_R_1_flag
    strange_PGCC["Dust_Temp_R_2_flag"]=Dust_R_2_flag
    strange_PGCC["CO_flag"]=CO_flag

    # Filter for new information
    strange_PGCC = strange_PGCC[(strange_PGCC.NHI_flag == 1) | (strange_PGCC.IRIS_3000_flag == 1) | (strange_PGCC.PLANCK_857_flag == 1)
                               | (strange_PGCC.PLANCK_545_flag == 1)| (strange_PGCC.PLANCK_353_flag == 1)| (strange_PGCC.Color_flag == 1)
                               | (strange_PGCC.Dust_Temp_R_1_flag == 1)| (strange_PGCC.Dust_Temp_R_2_flag == 1)
                               | (strange_PGCC.CO_flag == 1)]
    # Find amount of deviating parameters
    strangeness = [0]*len(strange_PGCC)
    for i in range(0, len(strange_PGCC)):
        strangeness[i] = strange_PGCC["NHI_flag"].iloc[i]+strange_PGCC["IRIS_3000_flag"].iloc[i]+strange_PGCC["PLANCK_857_flag"].iloc[i]+strange_PGCC["PLANCK_545_flag"].iloc[i]+strange_PGCC["PLANCK_353_flag"].iloc[i]+strange_PGCC["Color_flag"].iloc[i]+strange_PGCC["Dust_Temp_R_1_flag"].iloc[i]+strange_PGCC["Dust_Temp_R_2_flag"].iloc[i]
    strange_PGCC["strangeness"]=strangeness
    print(len(strange_PGCC),"/",len(pgcc), "deviate in at least one area")
    print(strangeness.count(1),"/",len(pgcc), "PGCC's deviate in one area")
    print(strangeness.count(2),"/",len(pgcc), "PGCC's deviate in two areas")
    print(strangeness.count(3),"/",len(pgcc), "PGCC's deviate in three areas")
    print(strangeness.count(4),"/",len(pgcc), "PGCC's deviate in four areas")
    print(strangeness.count(5),"/",len(pgcc), "PGCC's deviate in five areas")
    print(strangeness.count(6),"/",len(pgcc), "PGCC's deviate in six areas")
    print(strangeness.count(7),"/",len(pgcc), "PGCC's deviate in seven areas")
    print(strangeness.count(8),"/",len(pgcc), "PGCC's deviate in eight areas")
    print(strangeness.count(9),"/",len(pgcc), "PGCC's deviate in nine areas\n")

    #create flags
    flag = ["NHI_flag", "IRIS_3000_flag", "PLANCK_857_flag", "PLANCK_545_flag", "PLANCK_353_flag", "Color_flag", "Dust_Temp_R_1_flag", "Dust_Temp_R_2_flag", "CO_flag"]
    header_name = ["NHI", "IRIS_3000", "PLANCK_857", "PLANCK_545", "PLANCK_353", "Color", "Dust_Temp_R_1", "Dust_Temp_R_2", "CO"]

    for i in range(0, len(flag)):
        print(np.sum(strange_PGCC[flag[i]].values),"/",len(pgcc), "deviate ",header_name[i])

    ###
    # Deeper strangeness analysis. Iterate over strangeness
    ###
    if(show_figure==True):
        for i in range(0, 9):
            tmp_pgcc = strange_PGCC[strange_PGCC[flag[i]]==1]

            fig, axes = plt.subplots(nrows=3, ncols=3)
            fig.set_figheight(10)
            fig.set_figwidth(15)
            #axes[-1,-1].axis('off')
            x = ["NHI_avrg","IRIS_3000_avrg","PLANCK_857_avrg","PLANCK_545_avrg","PLANCK_353_avrg", "Color_avrg", "Dust_Temp_R_1_avrg","Dust_Temp_R_2_avrg", "CO_avrg"]*8
            y = [["Latitude", "IRIS_3000_avrg","PLANCK_857_avrg","PLANCK_545_avrg","PLANCK_353_avrg", "Color_avrg", "Dust_Temp_R_1_avrg","Dust_Temp_R_2_avrg", "CO_avrg"],
                 ["Latitude", "NHI_avrg","PLANCK_857_avrg","PLANCK_545_avrg","PLANCK_353_avrg", "Color_avrg", "Dust_Temp_R_1_avrg","Dust_Temp_R_2_avrg", "CO_avrg"],
                 ["Latitude", "NHI_avrg","IRIS_3000_avrg","PLANCK_545_avrg","PLANCK_353_avrg", "Color_avrg", "Dust_Temp_R_1_avrg","Dust_Temp_R_2_avrg", "CO_avrg"],
                 ["Latitude", "NHI_avrg","IRIS_3000_avrg","PLANCK_857_avrg","PLANCK_353_avrg", "Color_avrg", "Dust_Temp_R_1_avrg","Dust_Temp_R_2_avrg", "CO_avrg"],
                 ["Latitude", "NHI_avrg","IRIS_3000_avrg","PLANCK_857_avrg","PLANCK_545_avrg", "Color_avrg", "Dust_Temp_R_1_avrg","Dust_Temp_R_2_avrg", "CO_avrg"],
                 ["Latitude", "NHI_avrg","IRIS_3000_avrg","PLANCK_857_avrg","PLANCK_545_avrg", "PLANCK_353_avrg", "Dust_Temp_R_1_avrg","Dust_Temp_R_2_avrg", "CO_avrg"],
                 ["Latitude", "NHI_avrg","IRIS_3000_avrg","PLANCK_857_avrg","PLANCK_545_avrg", "PLANCK_353_avrg", "Color_avrg","Dust_Temp_R_2_avrg", "CO_avrg"],
                 ["Latitude", "NHI_avrg","IRIS_3000_avrg","PLANCK_857_avrg","PLANCK_545_avrg", "PLANCK_353_avrg", "Color_avrg","Dust_Temp_R_1_avrg", "CO_avrg"],
                 ["Latitude", "NHI_avrg","IRIS_3000_avrg","PLANCK_857_avrg","PLANCK_545_avrg", "PLANCK_353_avrg", "Color_avrg","Dust_Temp_R_1_avrg","Dust_Temp_R_2_avrg"]]
            x_names = ["N$_{HI}$ / cm$^{-2}$","$I_{3THz}$ / MJy/sr","$I_{857GHz}$ / MJy/sr","$I_{545GHz}$ / MJy/sr","$I_{353GHz}$ / MJy/sr", "Color $C$", "$T_{Dust,R1}$ / K","$T_{Dust,R2}$ / K", "W$_{CO}$ / K km/s"]*8
            y_names = [["Latitude / Degree", "$I_{3THz}$ / MJy/sr","$I_{857GHz}$ / MJy/sr","$I_{545GHz}$ / MJy/sr","$I_{353GHz}$ / MJy/sr", "Color $C$", "$T_{Dust,R1}$ / K","$T_{Dust,R2}$ / K", "W$_{CO}$ / K km/s"],
                 ["Latitude / Degree", "N$_{HI}$ / cm$^{-2}$","$I_{857GHz}$ / MJy/sr","$I_{545GHz}$ / MJy/sr","$I_{353GHz}$ / MJy/sr", "Color $C$", "$T_{Dust,R1}$ / K","$T_{Dust,R2}$ / K", "W$_{CO}$ / K km/s"],
                 ["Latitude / Degree", "N$_{HI}$ / cm$^{-2}$","$I_{3THz}$ / MJy/sr","$I_{545GHz}$ / MJy/sr","$I_{353GHz}$ / MJy/sr", "Color $C$", "$T_{Dust,R1}$ / K","$T_{Dust,R2}$ / K", "W$_{CO}$ / K km/s"],
                 ["Latitude / Degree", "N$_{HI}$ / cm$^{-2}$","$I_{3THz}$ / MJy/sr","$I_{857GHz}$ / MJy/sr","$I_{353GHz}$ / MJy/sr", "Color $C$", "$T_{Dust,R1}$ / K","$T_{Dust,R2}$ / K", "W$_{CO}$ / K km/s"],
                 ["Latitude / Degree", "N$_{HI}$ / cm$^{-2}$","$I_{3THz}$ / MJy/sr","$I_{857GHz}$ / MJy/sr","$I_{545GHz}$ / MJy/sr", "Color $C$", "$T_{Dust,R1}$ / K","$T_{Dust,R2}$ / K", "W$_{CO}$ / K km/s"],
                 ["Latitude / Degree", "N$_{HI}$ / cm$^{-2}$","$I_{3THz}$ / MJy/sr","$I_{857GHz}$ / MJy/sr","$I_{545GHz}$ / MJy/sr", "$I_{353GHz}$ / MJy/sr", "$T_{Dust,R1}$ / K","$T_{Dust,R2}$ / K", "W$_{CO}$ / K km/s"],
                 ["Latitude / Degree", "N$_{HI}$ / cm$^{-2}$","$I_{3THz}$ / MJy/sr","$I_{857GHz}$ / MJy/sr","$I_{545GHz}$ / MJy/sr", "$I_{353GHz}$ / MJy/sr", "Color $C$","$T_{Dust,R2}$ / K", "W$_{CO}$ / K km/s"],
                 ["Latitude / Degree", "N$_{HI}$ / cm$^{-2}$","$I_{3THz}$ / MJy/sr","$I_{857GHz}$ / MJy/sr","$I_{545GHz}$ / MJy/sr", "$I_{353GHz}$ / MJy/sr", "Color $C$","$T_{Dust,R1}$ / K", "W$_{CO}$ / K km/s"],
                 ["Latitude / Degree", "N$_{HI}$ / cm$^{-2}$","$I_{3THz}$ / MJy/sr","$I_{857GHz}$ / MJy/sr","$I_{545GHz}$ / MJy/sr", "$I_{353GHz}$ / MJy/sr", "Color $C$","$T_{Dust,R1}$ / K","$T_{Dust,R2}$ / K"]]
            print(x[i])
            head = pgcc.columns.values
            z=0
            for n in range(0,3):
                for m in range(0,3):
                    axes[n,m].scatter(x=pgcc[x[i]].values, y=pgcc[y[i][z]].values, label = "S: 0")
                    axes[n,m].scatter(x=tmp_pgcc[x[i]].values, y=tmp_pgcc[y[i][z]].values)
                    for s in range(1,10):
                        strangeness_pgcc = tmp_pgcc[tmp_pgcc["strangeness"]==s]
                        axes[n,m].scatter(x=strangeness_pgcc[x[i]].values, y=strangeness_pgcc[y[i][z]].values, label = "S: %d"%(int(s)))
                    axes[n,m].legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=False)
                    axes[n,m].set_xlabel(x_names[i])
                    axes[n,m].set_ylabel(y_names[i][z])
                    z += 1
                    if z == 9:
                        break
            fig.tight_layout()
            plt.savefig(filepath+name+"_"+"pgcc_"+x[i]+".png")
            plt.show()
    plt.close()
    return strange_PGCC, mean, std

def strangeness_quantifier(data):
    if len(data.loc[data["strangeness"]>=1])!=0:
        strange_data=data.loc[data["strangeness"]>=1]
        flag = ["NHI_flag", "IRIS_3000_flag", "PLANCK_857_flag", "PLANCK_545_flag","PLANCK_353_flag", "Color_flag", "Dust_Temp_R_1_flag", "Dust_Temp_R_2_flag", "CO_flag"]
        for i in range(0, len(flag)):
            current_flag_data=strange_data.loc[strange_data[flag[i]]==1]
            current_flag_data_length=len(current_flag_data)
            print("\n---\n"+flag[i]+" "+str(current_flag_data_length)+"\n")
            for j in range(0, len(flag)):
                if(i!=j):
                    current_length = len(current_flag_data.loc[current_flag_data[flag[j]]==1])
                    print(flag[j]+" "+str(current_length)+"("+str(np.round((current_length/current_flag_data_length)*100,2))+"%)\n")

def PGCC_CO_analysis(PGCC_info, strange_PGCC, selected_data):
    global filtepath
    print("\nPGCC_CO_analysis\n")

    #Filter nach CO
    CO_Schwelle = 1
    PGCC_CO = PGCC_info.loc[PGCC_info["CO_avrg"]>=CO_Schwelle]
    PGCC_CO_strange = strange_PGCC.loc[strange_PGCC["CO_avrg"]>=CO_Schwelle]
    PGCC_none_CO = PGCC_info.loc[PGCC_info["CO_avrg"]<=CO_Schwelle]
    PGCC_none_CO_strange = strange_PGCC.loc[strange_PGCC["CO_avrg"]<=CO_Schwelle]
    plt.figure()
    plt.scatter(PGCC_info["Dust_Temp_R_1_avrg"],PGCC_info["CO_avrg"])
    plt.scatter(PGCC_CO["Dust_Temp_R_1_avrg"],PGCC_CO["CO_avrg"], label="%d/%d PGCC's beinhalten CO"%(len(PGCC_CO), len(PGCC_info)))
    plt.grid()
    mean, std = PGCC_CO["Dust_Temp_R_1_avrg"].agg(['mean', 'std']).round(decimals=10)
    mean_CO, std_CO = PGCC_CO["CO_avrg"].agg(['mean', 'std']).round(decimals=10)
    if len(PGCC_CO)>100:
        plt.axvline(x=np.amin(PGCC_CO[PGCC_CO["CO_avrg"]>mean_CO]["Dust_Temp_R_1_avrg"]), color="black", linestyle = "dotted", label="min. Temp: %.1fK\nmax. Temp: %.1fK"%(np.amin(PGCC_CO[PGCC_CO["CO_avrg"]>mean_CO]["Dust_Temp_R_1_avrg"]),np.amax(PGCC_CO[PGCC_CO["CO_avrg"]>mean_CO]["Dust_Temp_R_1_avrg"])))
        plt.axvline(x=np.amax(PGCC_CO[PGCC_CO["CO_avrg"]>mean_CO]["Dust_Temp_R_1_avrg"]), color="black", linestyle = "dotted")
    else:
        plt.axvline(x=np.amin(PGCC_CO["Dust_Temp_R_1_avrg"]), color="black", linestyle = "dotted", label="min. Temp: %.3f K\nmax. Temp: %.3f K"%(np.amin(PGCC_CO["Dust_Temp_R_1_avrg"]),np.amax(PGCC_CO["Dust_Temp_R_1_avrg"])))
        plt.axvline(x=np.amax(PGCC_CO["Dust_Temp_R_1_avrg"]), color="black", linestyle = "dotted")
    plt.legend()
    plt.xlabel("$T_{Dust, R1}$ / K")
    plt.ylabel("$W_{CO}$ / K$\cdot$km/s")
    plt.savefig(filepath+"/CO_Dust_Temp_PGCC_Abhängigkeit"+name+filetype,bbox_inches='tight')
    plt.show()

    if(len(PGCC_CO)!=0):
        # Stelle die PGCC's dar


        print("\n%d / %d of the PGCC's contain CO\n%d / %d of them are strange'\n\n"%(int(len(PGCC_CO)),int(len(PGCC_info)), int(len(PGCC_CO_strange)),int(len(PGCC_info))))

        # Bestimme die Durchschnittswerte der PGCC Temperatur die CO Spuren aufweisen
        mean = [0]*4
        std = [0]*4
        mean[0], std[0] = PGCC_CO["Dust_Temp_R_1_avrg"].agg(['mean', 'std']).round(decimals=10)
        mean[1], std[1] = PGCC_CO_strange["Dust_Temp_R_1_avrg"].agg(['mean', 'std']).round(decimals=10)
        mean[2], std[2] = PGCC_CO["Dust_Temp_R_2_avrg"].agg(['mean', 'std']).round(decimals=10)
        mean[3], std[3] = PGCC_CO_strange["Dust_Temp_R_2_avrg"].agg(['mean', 'std']).round(decimals=10)
        print("Dust_Temp_R1 (all PGCC's):\nAverage temperature %.6f, standard deviation %.6f\n"%(mean[0], std[0]))
        print("Dust_Temp_R2 (all PGCC's):\nAverage temperature %.6f, standard deviation %.6f\n"%(mean[2], std[2]))
        print("Duste_Temp_R1 (strange PGCC's):\nAverage temperature %.6f, standard deviation %.6f\n"%(mean[1], std[1]))
        print("Dust_Temp_R2 (strange PGCC's):\nAverage temperature %.6f, standard deviation %.6f\n"%(mean[3], std[3]))

        #Intensitäte im FIR und NHI darstellen
        fig, axes = plt.subplots(nrows=1, ncols=2)
        fig.set_figheight(5)
        fig.set_figwidth(10)
        mean, std= PGCC_CO["IRIS_3000_avrg"].agg(['mean', 'std']).round(decimals=10)
        axes[0].hist(PGCC_CO["IRIS_3000_avrg"], bins=181, log=True, label ="$\mu=%.3f\n \sigma=%.3f$"%(mean, std))
        axes[0].set_xlabel("$I_{3THz}$")
        mean, std= PGCC_CO["NHI_avrg"].agg(['mean', 'std']).round(decimals=10)
        axes[1].hist(PGCC_CO["NHI_avrg"], bins=181, log=True, label ="$\mu=%.3f\n \sigma=%.3f$"%(mean, std))
        axes[1].set_xlabel("$NHI$")
        plt.show()

        # Find CO-Dark, -Bright and -Strong PGCC's and visualize them
        m,b,c=korrelation(PGCC_info, "NHI_avrg", "CO_avrg", "lin_fit", "test", "test", False)
        pgcc_dark = PGCC_info.copy()
        drop_list=[]
        for i in range(0,len(PGCC_info)):
            if(PGCC_info["CO_avrg"].iloc[i]<m*PGCC_info["NHI_avrg"].iloc[i]+b+3*std):
                drop_list.append(i)

        pgcc_super_bright = pgcc_dark.drop(pgcc_dark.index[drop_list])
        pgcc_bright = pd.concat([PGCC_info[PGCC_info["CO_avrg"]>1], pgcc_super_bright,  pgcc_super_bright]).drop_duplicates(keep=False)
        pgcc_dark = PGCC_info[PGCC_info["CO_avrg"]<=1]
        plt.scatter(pgcc_bright["NHI_avrg"],pgcc_bright["CO_avrg"],c="steelblue", label="%d - Bright PGCC's"%(len(pgcc_bright)))
        plt.scatter(pgcc_super_bright["NHI_avrg"],pgcc_super_bright["CO_avrg"],c="orange", label="%d - CO-strong PGCC's"%(len(pgcc_super_bright)))
        plt.scatter(pgcc_dark["NHI_avrg"],pgcc_dark["CO_avrg"],c="black", label="%d - Dark PGCC"%(len(pgcc_dark)))#pgcc_dark["Dust_Temp_R_1_avrg"])
        plt.xlabel(global_labels[0])
        plt.ylabel(global_labels[-1])
        plt.grid()
        plt.legend()
        plt.savefig(filepath+"/CO_Dark_Bright"+name+filetype,bbox_inches='tight')
        plt.show()

        plt.figure(figsize=(21, 12))
        wideview(pgcc_dark["Longitude"], pgcc_dark["Latitude"], "black", 1, "%d/%d - Dark PGCC's"%(len(pgcc_dark),len(PGCC_info)))
        wideview(pgcc_bright["Longitude"], pgcc_bright["Latitude"], "steelblue", 1, "%d/%d - Bright PGCC's"%(len(pgcc_bright),len(PGCC_info)))
        wideview(pgcc_super_bright["Longitude"], pgcc_super_bright["Latitude"], "orange", 1, "%d/%d - CO-leuchtkräfigte PGCC's"%(len(pgcc_super_bright),len(PGCC_info)))
        plt.legend()
        plt.xlabel("Longitude / grad")
        plt.ylabel("Latitude / grad")
        plt.savefig(filepath+"/wideview_CO_Bright_dark"+filetype,bbox_inches='tight')
        plt.show()



    else:
        print("The chosen are doesn't contains CO emitting PGCC's")

def TWO_D_HISTO(x,y,z, xlabel, ylabel, zlabel):
    plt.scatter(x, y, c=z)
    x_fig = xlabel
    y_fig = ylabel
    for i in range(0, len(global_labels)):
        if xlabel in global_labels[i]:
            xlabel=global_labels[i]
        if ylabel in global_labels[i]:
            ylabel=global_labels[i]
        if zlabel in global_labels[i]:
            zlabel=global_labels[i]
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()
    plt.colorbar(label=zlabel)
    plt.savefig(filepath+name+"_"+x_fig+"_"+y_fig+"_histo.png")
    plt.show()


def Vergleiche(data1, data2, name1, name2):
    # PGCC's
    print("PGCC-Comparison\n\n%s: %d PGCC's \n%s: %d PGCC's"%(name1,len(data1),name2,len(data2)))
    header = ["NHI_avrg", "IRIS_3000_avrg", "PLANCK_857_avrg", "PLANCK_545_avrg", "PLANCK_353_avrg", "Color_avrg", "Dust_Temp_R_1_avrg", "Dust_Temp_R_2_avrg", "CO_avrg"]
    mean1 = [0]*9 #NHI - IRIS - 3xPLANCK - Color - Temp1 -Temp2 -CO
    std1 = [0]*9
    mean2 = [0]*9
    std2 = [0]*9
    for i in range(0, len(header)):
        mean1[i], std1[i] = data1[header[i]].agg(['mean', 'std']).round(decimals=10)
        mean2[i], std2[i] = data2[header[i]].agg(['mean', 'std']).round(decimals=10)
        print("\n%s"%(header[i]))
        print("%s: mean %.5f, std %.5f\n%s: mean %.5f, std %.5f"%(name1, mean1[i],std1[i],name2, mean2[i],std2[i]))
    CO_Schwelle = 0.5
    print("\n%s: PGCC's mit CO gehalt - %d / %d\n%s: PGCC's mit CO gehalt - %d / %d"%(name1, len(data1.loc[data1["CO_avrg"]>=CO_Schwelle]), len(data1),name2, len(data2.loc[data2["CO_avrg"]>=CO_Schwelle]), len(data2)))
    x_label = ["Latitude / Degree", "NHI / cm$^{-2}$", "NHI / cm$^{-2}$", "NHI / cm$^{-2}$", "NHI / cm$^{-2}$", "NHI / cm$^{-2}$", "Latitude / Degree", "Latitude / Degree", "Latitude / Degree"]
    y_label = ["NHI / cm$^{-2}$",  "$I_{3THz}$ / MJycm$^2$/sr", "$I_{857GHz}$ / MJycm$^2$/sr", "$I_{545GHz}$ / MJycm$^2$/sr", "$I_{353GHz}$ / MJycm$^2$/sr", "Color $C$", "$T_{Dust, R1}$ / K", "$T_{Dust, R2}$ / K", "$W_{CO}$ / K $\cdot$ km/s"]
    x_name = ["Latitude", "NHI_avrg", "NHI_avrg", "NHI_avrg", "NHI_avrg", "NHI_avrg", "Latitude", "Latitude", "Latitude"]
    y_name = ["NHI_avrg",  "IRIS_3000_avrg", "PLANCK_857_avrg", "PLANCK_545_avrg", "PLANCK_353_avrg", "Color_avrg", "Dust_Temp_R_1_avrg", "Dust_Temp_R_2_avrg", "CO_avrg"]
    z = 0
    fig, axes = plt.subplots(nrows=3, ncols=3)
    fig.set_figheight(10)
    fig.set_figwidth(15)
    for i in range(0,3):
        for j in range(0,3):
            axes[i,j].scatter(x=data1[x_name[z]],y=data1[y_name[z]], label="%s"%(name1))
            axes[i,j].scatter(x=data2[x_name[z]],y=data2[y_name[z]], label="%s"%(name2))
            axes[i,j].set_xlabel(x_label[z])
            axes[i,j].set_ylabel(y_label[z])
            axes[i,j].legend(loc="upper right")
            z+=1
    fig.tight_layout()
    plt.savefig(filepath+"/Comparison_"+name1+"_"+name2+".png")
    plt.show()

    fig, axes = plt.subplots(nrows=3, ncols=3)
    fig.set_figheight(10)
    fig.set_figwidth(15)
    #axes[-1,-1].axis('off')
    mean = [0]*9
    std = [0]*9
    z = 2
    head = data1.columns.values
    for i in range(0,3):
        for j in range(0,3):
            mean[z-2], std[z-2] = data1[head[z]].agg(['mean', 'std']).round(decimals=10)
            data1[head[z]].plot.hist(ax=axes[i,j],grid=True, bins=100, label ="%s\n$\mu$=%.6f\n$\sigma$=%.6f"%(name1,mean[z-2], std[z-2]),subplots=True)
            data2[head[z]].plot.hist(ax=axes[i,j],grid=True, bins=100, label =name2,subplots=True, color="orange")
            axes[i,j].legend()
            axes[i,j].set_xlabel(global_labels_no_pgcc[z-2])
            axes[i,j].set_title(global_names[z-2])
            plt.grid(axis="x", alpha=0.75)
            #data.plot.scatter(ax=axes[i,j],x="Longitude",y="Latitude",c=head[z],s=1,subplots=True, colormap='viridis')
            z += 1
            if z == 11: break;
    fig.tight_layout()
    plt.savefig(filepath+name1+"_"+name2+"_histo_vgl.png")
    plt.show()


def PGCC_slice_analysator(data, strange_data):
    slices = len(data)
    total_pgcc_count = 0
    for i in range(0,slices):
        total_pgcc_count += len(data[i])
    x = []
    x_strange = []
    y = []
    y_strange = []
    c = []
    c_strange = []
    for i in range(0, slices):
        if len(data[i]!=0):
            x.append(len(data[i]))
            y.append(len(data[i])/total_pgcc_count)
            c.append((np.amax(data[i]["Latitude"])+np.amin(data[i]["Latitude"]))/2)
            x_strange.append(len(strange_data[i]))
            y_strange.append(len(strange_data[i])/len(data[i]))
            c_strange.append((np.amax(strange_data[i]["Latitude"])+np.amin(strange_data[i]["Latitude"]))/2)
        else:
            x.append(0)
            y.append(0)
            c.append(0)
            x_strange.append(0)
            y_strange.append(0)
            c_strange.append(0)
    plt.figure()
    plt.scatter(x,y,c=c)
    plt.colorbar(label="Latitude")
    plt.xlabel("PGCC's found")
    plt.ylabel("Share of PGCC's in the region")
    plt.show()
    plt.scatter(x_strange,y_strange,c=c_strange)
    plt.colorbar(label="Latitude")
    plt.xlabel("strange PGCC's found")
    plt.ylabel("Share of strange PGCC's in the region")
    plt.show()

####
# Function to create two overlaying plots with different y-axis without errorbars
####
def double_plot(x,xlabel,y1,y1_label,y2,y2_label,alpha,invert_y1, invert_y2,filename):
    fig, ax1 = plt.subplots()
    color = 'tab:orange'
    ax1.set_xlabel('Latitude / $\circ$')
    ax1.set_ylabel(y1_label, color="black")
    if invert_y1 == True:
        ax1.invert_yaxis()
    ax1.scatter(x, y1, color=color, alpha=alpha)
    ax1.tick_params(axis='y', labelcolor=color)
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:blue'
    ax2.set_ylabel(y2_label, color="black")  # we already handled the x-label with ax1
    ax2.scatter(x, y2, color=color, alpha=alpha)
    ax2.tick_params(axis='y', labelcolor=color)
    if invert_y2 == True:
        ax2.invert_yaxis()
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.grid()
    plt.savefig(filepath+filename+filetype,bbox_inches='tight')
    plt.show()

####
# Function to create two overlaying plots with different y-axis and errorbars
####
def double_plot_yerr(x,xlabel,y1,y1err,y1_label,y2,y2err,y2_label,alpha,invert_y1, invert_y2,filename):
    fig, ax1 = plt.subplots()
    color = 'tab:orange'
    ax1.set_xlabel('Latitude / $\circ$')
    ax1.set_ylabel(y1_label, color="black")
    ax1.set_ylim(np.amin(y1)-0.5,np.amax(y1)+0.5)
    if invert_y1 == True:
        ax1.invert_yaxis()
    print(x)
    ax1.scatter(x=x,y=y1,alpha=alpha,color=color)
    ax1.errorbar(x=x, y=y1,yerr=y1err, color=color, alpha=alpha, linestyle="None")
    ax1.tick_params(axis='y', labelcolor=color)
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel(y2_label, color="black")
    ax2.set_ylim(np.amin(y2)-0.05,np.amax(y2)+0.05)
    ax2.scatter(x,y2,color=color)
    ax2.errorbar(x=x, y=y2,yerr=y2err, color=color, alpha=alpha, linestyle="None")
    ax2.tick_params(axis='y', labelcolor=color)
    if invert_y2 == True:
        ax2.invert_yaxis()
    fig.tight_layout()
    plt.grid()
    plt.savefig(filepath+filename+"_yerr"+filetype,bbox_inches='tight')
    plt.show()


####
# Analysis tool for a skyselection - It includes a general sky-, PGCC-, neighbourhood- and CO-analysis
####
def sky_analyser(all_sky_table,name_input,lon_lim_min, lon_lim_max, lat_lim_min, lat_lim_max, inverted, lon_steps, lat_steps, filter_option, show_figure, debugging, dist):
    global filepath
    global name
    name = name_input
    filepath = "./output/"
    data_slices = []
    if os.path.exists(filepath+name):
        shutil.rmtree(filepath+name)
    os.makedirs(filepath+name)
    filepath=filepath+name+"/"
    print("Your files will be saved in "+filepath)
    all_sky_table=all_sky_table
    selected_data, data_slices = data_selector(all_sky_table,lon_lim_min, lon_lim_max, lat_lim_min, lat_lim_max, inverted, lon_steps, lat_steps, debugging)

    print("\n...Dataset generated...\n")
    if filter_option == False:
        test_map = debugging_tools.gen_test_map(selected_data, "NHI",debugging, 0)
    if filter_option != False:
        selected_data=filter_data(selected_data, filter_option)
        if(lat_steps!= 0 or lon_steps != 0):
            for i in range(0, lat_steps):
                data_slices[i]=filter_data(data_slices[i], filter_option)
        print("\n...Dataset filtered...\n")
    if((show_figure == True) and (np.abs(lon_lim_max-lon_lim_min)<=100) and (np.abs(lat_lim_max)-np.abs(lat_lim_min)<=100)):
        show_pics(data_slices,lon_steps, lat_steps, False, False)
    m_complete,b_complete,c_complete=korrelation(selected_data, "NHI", "IRIS_3000GHz", "lin_fit", "NHI / cm$^-2$", "$I_{FIR}$ / MJy/sr", show_figure)
    m = []
    b = []
    c = []
    if(lat_steps != 0):
        for i in range(0, len(data_slices)):
            m.append([])
            b.append([])
            c.append([])
            m[i],b[i],c[i]=korrelation(data_slices[i], "NHI", "IRIS_3000GHz", "lin_fit", "NHI", "$I_{FIR}$", False)
    NHI_avrg = mitteln(data_slices, "NHI", lon_steps, lat_steps)
    TMP_R1_avrg = mitteln(data_slices, "Dust_Temp_R1", lon_steps, lat_steps)
    TMP_R2_avrg = mitteln(data_slices, "Dust_Temp_R2", lon_steps, lat_steps)
    Color_avrg = mitteln(data_slices, "Color", lon_steps, lat_steps)
    Latitude_avrg = mitteln(data_slices, "Latitude", lon_steps, lat_steps)
    if(show_figure==True and lat_steps != 0):
        print("---!!!The following plots contain the complete ISM data and not only PGCC's!!!---")
        scatter_plot(Latitude_avrg, NHI_avrg, "Latitude / Degree", "N$_{HI}$ / cm$^{-2}$","ISM_Latitude_NHI", lon_steps, lat_steps, True)
        scatter_plot(Latitude_avrg, c, "Latitude / Degree", "Correlationcoeff. of N$_{HI}$-$I_{3THz}$","ISM_correlation_of_NHI_I3000GHz", lon_steps, lat_steps, True)
        scatter_plot(Latitude_avrg, TMP_R1_avrg, "Latitude / Degree", "$T_{Dust, R1}$ / K","ISM_Laitude_Dust_Temp_R1", lon_steps, lat_steps, True)
        scatter_plot(Latitude_avrg, TMP_R2_avrg, "Latitude / Degree", "$T_{Dust,R2}$ / K","ISM_Latitude_Dust_Temp_R2", lon_steps, lat_steps, True)
        double_plot(Latitude_avrg,"Latitude / Degree",TMP_R1_avrg,"$T_{Dust,R1}$ / K",Color_avrg,"ISM_Color C",0.4,False,True,"doubleplot_Color_TempR1_"+name)
        double_plot(Latitude_avrg,"Latitude / Degree",TMP_R2_avrg,"$T_{Dust,R2}$ / K",Color_avrg,"ISM_Color C",0.4,False,True,"doubleplot_Color_TempR2_"+name)
        fig, axes = plt.subplots(nrows=1, ncols=2)
        fig.set_figheight(5)
        fig.set_figwidth(10)
        mean_R1 = []
        std_R1 = []
        mean_R2 = []
        std_R2 = []
        for i in range(0,len(data_slices)):
            mean, std = data_slices[i]["Dust_Temp_R1"].agg(['mean','std']).round(decimals=10)
            mean_R1.append(mean)
            std_R1.append(std)
            if(np.sum(data_slices[i]["Latitude"])>=0):
                data_slices[i]["Dust_Temp_R1"].plot.hist(ax=axes[0],grid=True, bins=100, label ="Lat:%.1f$^\circ$ - %.1f$^\circ$ \nMean=%.1fK\nSTD=%.1fK"%(np.round(np.amin(data_slices[i]["Latitude"]),0),np.round(np.amax(data_slices[i]["Latitude"]),0),mean, std))
                axes[0].set_xlabel("Dusttemperature $T_{Dust, R1}$")
            else:
                data_slices[i]["Dust_Temp_R1"].plot.hist(ax=axes[1],grid=True, bins=100, label ="Lat:%.1f$^\circ$ - %.1f$^\circ$ \nMean=%.1fK\nSTD=%.1fK"%(np.round(np.amin(data_slices[i]["Latitude"]),0),np.round(np.amax(data_slices[i]["Latitude"]),0),mean, std),zorder=-i)
                axes[1].set_xlabel("Dusttemperature $T_{Dust, R1}$")
        axes[0].legend(loc="upper right")
        axes[1].legend(loc="upper right")
        plt.savefig(filepath+name+"ISM_Latitude_Temp_R1.png")
        plt.show()
        fig, axes = plt.subplots(nrows=1, ncols=2)
        fig.set_figheight(5)
        fig.set_figwidth(10)
        for i in range(0,len(data_slices)):
            mean, std = data_slices[i]["Dust_Temp_R2"].agg(['mean','std']).round(decimals=10)
            mean_R2.append(mean)
            std_R2.append(std)
            if(np.sum(data_slices[i]["Latitude"])>=0):
                data_slices[i]["Dust_Temp_R2"].plot.hist(ax=axes[0],grid=True, bins=100, label ="Lat:%.1f$^\circ$ - %.1f$^\circ$ \nMean=%.1fK\nSTD=%.1fK"%(np.round(np.amin(data_slices[i]["Latitude"]),0),np.round(np.amax(data_slices[i]["Latitude"]),0),mean, std))
                axes[0].set_xlabel("Dusttemperature $T_{Dust, R1}$")
            else:
                data_slices[i]["Dust_Temp_R2"].plot.hist(ax=axes[1],grid=True, bins=100, label ="Lat:%.1f$^\circ$ - %.1f$^\circ$ \nMean=%.1fK\nSTD=%.1fK"%(np.round(np.amin(data_slices[i]["Latitude"]),0),np.round(np.amax(data_slices[i]["Latitude"]),0),mean, std),zorder=-i)
                axes[1].set_xlabel("Dusttemperature $T_{Dust, R1}$")
        axes[0].legend(loc="upper right")
        axes[1].legend(loc="upper right")
        plt.savefig(filepath+name+"ISM_Latitude_Temp_R2.png")
        plt.show()
        plt.errorbar(x=Latitude_avrg, y=mean_R1, yerr=std_R1, linestyle="None")
        mean = pd.DataFrame(data={"Latitude":Latitude_avrg,"mean_R1":mean_R1})
        slope, intercept, r_value, p_value, std_err = stats.linregress(mean[mean["Latitude"]>0]["Latitude"], mean[mean["Latitude"]>0]["mean_R1"])
        popt,perr=anpassung_yerr(func_lin, mean[mean["Latitude"]>0]["Latitude"],mean[mean["Latitude"]>0]["mean_R1"],std_R1, [0,0], False, False)
        plt.plot(np.linspace(np.amin(mean[mean["Latitude"]>0]["Latitude"]), np.amax(mean[mean["Latitude"]>0]["Latitude"]),100), popt[0]*np.linspace(np.amin(mean[mean["Latitude"]>0]["Latitude"]), np.amax(mean[mean["Latitude"]>0]["Latitude"]),100)+popt[1], label = "T(h)=$(%.3f\pm %.3f) K/^\circ \cdot h+(%.3f\pm %.3f) K$\n$std err=%.3f$"%(popt[0],perr[0],popt[1],perr[1],std_err), color="orange")
        plt.grid()
        plt.legend()
        plt.xlabel("Latitude")
        plt.ylabel("Dusttemperature $T_{Dust, R1}$")
        plt.savefig(filepath+"/ISM_Temp_Breitenabhängigkeit_"+name+".png",bbox_inches='tight')
        plt.show()
        plt.errorbar(x=Latitude_avrg, y=np.log(mean_R1), yerr=std_R1, linestyle="None")
        plt.xlabel("Latitude")
        plt.ylabel("log $T_{Dust, R1}$")
        plt.show()
        scatter_plot(Latitude_avrg, Color_avrg, "Latitude / Degree", "Color","ISM_Latitude_color", lon_steps, lat_steps, True)

    ####
    # Begin PGCC-Analysis
    ####
    print("\n---PGCC-Analysis---\n")
    PGCC_selected = selected_data.loc[selected_data['PGCC'] == 1]
    PGCC_info = PGCC_eight_pix(all_sky_table,PGCC_selected["Longitude"], PGCC_selected["Latitude"])
    if(show_figure == True):
        scatter_plot(PGCC_info["Dust_Temp_R_1_avrg"].values, PGCC_info["NHI_avrg"].values, "$T_{Dust, R1}$ / K","N$_{HI}$ / cm$^{-2}$","Dust_Temp_R1_NHI_PGCC", lon_steps, lat_steps, True)
        scatter_plot(PGCC_info["NHI_avrg"].values, PGCC_info["Color_avrg"].values, "N$_{HI}$ / cm$^{-2}$", "Color","NHI_Color_PGCC", lon_steps, lat_steps, True)
        scatter_plot(PGCC_info["NHI_avrg"].values, PGCC_info["Dust_Temp_R_1_avrg"].values, "N$_{HI}$ / cm$^{-2}$", "$T_{Dust, R1}$ / K","NHI_Dust_Temp_R1_PGCC", lon_steps, lat_steps, True)
        scatter_plot(PGCC_info["Latitude"].values, PGCC_info["NHI_avrg"].values, "Latitude / Degree", "N$_{HI}$ / cm$^{-2}$","Latitude_NHI_PGCC", lon_steps, lat_steps, True)
        scatter_plot(PGCC_info["Latitude"].values, PGCC_info["Color_avrg"].values, "Latitude / Degree", "Color","Latitude_Color_PGCC", lon_steps, lat_steps, True)
        scatter_plot(PGCC_info["Latitude"].values, PGCC_info["Dust_Temp_R_1_avrg"].values, "Latitude", "$T_{Dust, R1}$ / K", "Latitude_Dust_Temp_R1_PGCC", lon_steps, lat_steps, True)
    ####
    # PGCC-Neighbourhood Analysis
    ####
    print("\n---Neighbourhood-analysis---\n")
    PGCC_pairs, PGCC_expanded = find_PGCCs_nearby(PGCC_info, dist,0.2)
    if (PGCC_pairs!=False):
        expanded_pgcc = [item for elem in PGCC_expanded for item in elem]
        expanded_pgcc = list(dict.fromkeys(expanded_pgcc))
        PGCC_objects_list = [item for elem in PGCC_pairs for item in elem]
        no_pgcc_list = (np.setdiff1d(PGCC_objects_list,expanded_pgcc)).tolist()
    if PGCC_pairs == False:
        no_pgcc_list = False
    print("Neighbourhood data collected")
    if ((lon_steps==0 or lat_steps==0) and (np.abs(lon_lim_max-lon_lim_min)<=100) and (np.abs(lat_lim_max)-np.abs(lat_lim_min)<=100)):
        show_pics(data_slices,lon_steps, lat_steps, PGCC_pairs, no_pgcc_list)

    ####
    # advanced PGCC- and CO-analysis
    ####
    if len(PGCC_info)>1:
        print("\n---All selected PGCCs---\n")
        strange_PGCC, mean, std=pgcc_analyser(PGCC_info, PGCC_info, show_figure)
        strangeness_quantifier(strange_PGCC)
        PGCC_CO_analysis(PGCC_info, strange_PGCC, selected_data)
    if(lat_steps != 0):
        PGCC_slices = []
        PGCC_slices_info = []
        strange_PGCC_slices=[]
        for i in range(0, int(len(data_slices))):
            PGCC_slices.append([])
            PGCC_slices_info.append([])
            strange_PGCC_slices.append([])
            PGCC_slices[i]=data_slices[i].loc[data_slices[i]["PGCC"]==1]
            PGCC_slices_info[i]=PGCC_eight_pix(all_sky_table,PGCC_slices[i]["Longitude"], PGCC_slices[i]["Latitude"])
            print("\nSLICE_INFO:\nLat_min %.3f\nLat_max %.3f"%(np.amin(PGCC_slices_info[i]["Latitude"]), np.amax(PGCC_slices_info[i]["Latitude"])))
            strange_PGCC_slices[i], mean, std=pgcc_analyser(PGCC_slices_info[i], PGCC_info, False)
        mean_R1 = []
        std_R1 = []
        mean_R2 = []
        std_R2 = []
        mean_color = []
        std_color = []
        for i in range(0,len(data_slices)):
            mean, std = PGCC_slices_info[i]["Dust_Temp_R_1_avrg"].agg(['mean','std']).round(decimals=10)
            mean_R1.append(mean)
            std_R1.append(std)
            mean, std = PGCC_slices_info[i]["Dust_Temp_R_2_avrg"].agg(['mean','std']).round(decimals=10)
            mean_R2.append(mean)
            std_R2.append(std)
            mean_c, std_c = PGCC_slices_info[i]["Color_avrg"].agg(['mean','std']).round(decimals=10)
            mean_color.append(mean_c)
            std_color.append(std_c)
        PGCC_slice_analysator(PGCC_slices_info,strange_PGCC_slices)
        PGCC_Latitude_avrg = mitteln(PGCC_slices_info, "Latitude", lon_steps, lat_steps)
        double_plot(PGCC_Latitude_avrg,"Latitude / Degree",mean_R1,"$T_{Dust,R1} / K$",mean_color,"Color C",0.4,False,True,"doubleplot_Color_TempR1_PGCC_"+name)
        double_plot(PGCC_Latitude_avrg,"Latitude / Degree",mean_R2,"$T_{Dust,R2} / K$",mean_color,"Color C",0.4,False,True,"doubleplot_Color_TempR2_PGCC_"+name)
        plt.errorbar(x=PGCC_Latitude_avrg, y=mean_R1, yerr=std_R1, linestyle="None")
        mean = pd.DataFrame(data={"Latitude":PGCC_Latitude_avrg,"mean_R1":mean_R1})
        slope, intercept, r_value, p_value, std_err = stats.linregress(mean[mean["Latitude"]>0]["Latitude"], mean[mean["Latitude"]>0]["mean_R1"])
        popt,perr=anpassung_yerr(func_lin, mean[mean["Latitude"]>0]["Latitude"],mean[mean["Latitude"]>0]["mean_R1"],std_R1, [0,0], False, False)
        plt.plot(np.linspace(np.amin(mean[mean["Latitude"]>0]["Latitude"]), np.amax(mean[mean["Latitude"]>0]["Latitude"]),100), popt[0]*np.linspace(np.amin(mean[mean["Latitude"]>0]["Latitude"]), np.amax(mean[mean["Latitude"]>0]["Latitude"]),100)+popt[1], label = "T(h)=$(%.3f\pm %.3f) K/^\circ \cdot h+(%.3f\pm %.3f) K$\n$std err=%.3f$"%(popt[0],perr[0],popt[1],perr[1],std_err), color="orange")
        plt.legend()
        plt.grid()
        plt.xlabel("Latitude / Degree")
        plt.ylabel("$T_{Dust, R1}$")
        plt.savefig(filepath+"/PGCC_Temp_distribution"+name+filetype,bbox_inches='tight')
        plt.show()
    selected_data.to_csv(filepath+name+".csv", sep=";")
    print("---Analysis done---")
    return PGCC_info, strange_PGCC, PGCC_pairs, PGCC_expanded, no_pgcc_list
