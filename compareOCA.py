#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Created on 2022-09-06

Copyright (c) 2022 Erik Johansson

@author:     Erik Johansson
@contact:    <erik.johansson@smhi.se>
 
'''

import os
import sys
import numpy as np
import glob
import time
import pdb
import netCDF4 as nc  # @UnresolvedImport
import h5py
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from pyresample import load_area, save_quicklook, SwathDefinition
from pyresample.kd_tree import resample_nearest


def getMSGNum(y, m):
    #: MSG1
    if y in [2004, 2005, 2006]:
        mn = 1
    #: MSG1 / MSG2
    elif y in [2007]:
        if m in [1, 2, 3]:
            mn = 1
        #: Split month
        elif m in [4]:
            mn = 12
        else:
            mn = 2
    #: MSG2
    elif y in [2008, 2009, 2010, 2011, 2012]:
        mn = 2
    #: MSG1 / MSG2 / MSG3
    elif y in [2013]:
        #: Split month
        if m in [1]:
            mn = 23
        else:
            mn = 3
    #: MSG3
    elif y in [2014, 2015, 2016, 2017]:
        mn = 3
    #: MSG3 / MSG4
    elif y in [2018]:
        if m in [1]:
            mn = 3
        #: Split month
        elif m in [2]:
            mn = 34
        else:
            mn = 4
    #: MSG4
    elif y in [2019]:
        if m <= 8: 
            mn = 4
    else:
        mn = 0
    return mn

# /W_XX-EUMETSAT-Darmstadt,OCA+,MET09+SEVIRI_C_EUMG_20070414110000_1_OR_FES_E0000_0100.nc
# 
#         0100.nc
# 2007/11/04/W_XX-EUMETSAT-Darmstadt,OCA+,MET09+SEVIRI_C_EUMG_20071104121500_1_OR_FES_E0000_0100.nc
# 2007/10/04/W_XX-EUMETSAT-Darmstadt,OCA+,MET09+SEVIRI_C_EUMG_20071004184500_1_OR_FES_E0000_0100.nc'
# 2007/09/04/W_XX-EUMETSAT-Darmstadt,OCA+,MET09+SEVIRI_C_EUMG_20070904153000_1_OR_FES_E0000_0100.nc'
# 2007/08/04/W_XX-EUMETSAT-Darmstadt,OCA+,MET09+SEVIRI_C_EUMG_20070804194500_1_OR_FES_E0000_0100.nc'
# 2007/07/04/W_XX-EUMETSAT-Darmstadt,OCA+,MET09+SEVIRI_C_EUMG_20070704171500_1_OR_FES_E0000_0100.nc'
# 2007/06/04/W_XX-EUMETSAT-Darmstadt,OCA+,MET09+SEVIRI_C_EUMG_20070604123000_1_OR_FES_E0000_0100.nc'
# 2007/05/04/W_XX-EUMETSAT-Darmstadt,OCA+,MET09+SEVIRI_C_EUMG_20070504153000_1_OR_FES_E0000_0100.nc'



def compareOCAListname(fn, flist):
    same = 0
    dx = os.path.basename(fn).split('_')[4]
    i = -1
    for f in flist:
        i = i + 1
        d = os.path.basename(f).split('_')[4]
        if dx == d:
            same = 1
            break
    return same, i
                    

def getMeanOca(od, ys, ms, cm=False, lt=True):
    #: od = dir
    #: ys = year
    #: ms = month
    #: cm = create monthly mean
    #: lt = load temp files

    tempname = os.path.join('TempFiles', 'OCAMean', 'oca_monthly-mean_%d-%02d.h5' %(ys, ms))
    corruptFiles = []
    #===========================================================================
    # #: Calc num of satellites and num of files
    # msgwd = ''
    # monthNames = []
    # reserve = []
    # for x in range(4, 0 ,-1):
    #     monthNames_x = glob.glob('%s/MSG%d/%d/%02d/*/*.nc' %(od, x, ys, ms))
    #     if len(monthNames_x) != 0:
    #         msgwd = msgwd + '%s' %x
    #         if len(monthNames) == 0:
    #             monthNames = monthNames_x
    #         else:
    #             for fx in monthNames_x:
    #                 same, same_p = compareOCAListname(fx, monthNames)  # @UnusedVariable
    #                 if same == 1:
    #                     reserve.append(fx)
    #                 else:
    #                     monthNames.append(fx)
    # msgwd = msgwd[::-1]
    # #: Rerun for month with more than one satellite
    # if (len(msgwd) > 1):
    #     print(msgwd)
    #     cm = True
    #===========================================================================
    #===========================================================================
    # #: Update h5file
    # if os.path.isfile(tempname):
    #     h5f = h5py.File(tempname, 'a')
    #     if 'MSG Name' not in h5f.keys():
    #         h5f.create_dataset('MSG Name', data = [int(msgwd)])
    #     else:
    #         del h5f['MSG Name']
    #         h5f.create_dataset('MSG Name', data = [int(msgwd)])
    #     if 'Number of Files' not in h5f.keys():
    #         h5f.create_dataset('Number of Files', data = [len(monthNames)])
    #     else:
    #         del h5f['Number of Files']
    #         h5f.create_dataset('Number of Files', data = [len(monthNames)])
    #     if 'Corrupt Files' not in h5f.keys():
    #         h5f.create_dataset('Corrupt Files', data = corruptFiles)
    #     if 'Extra Files' not in h5f.keys():
    #         h5f.create_dataset('Extra Files', data = reserve)
    #     h5f.close()
    #===========================================================================
        
    if (cm == True) or (lt == False) or (not os.path.isfile(tempname)):
        print('Calculate OCA monthly mean')
        #: Get file names
        if False:
            monthNames_n = glob.glob('%s/*/*000000_1_OR*.nc' %od)
            monthNames_d = glob.glob('%s/*/*120000_1_OR*.nc' %od)
            monthNames = []
            for d in [5, 10]:
                monthNames.append(monthNames_n[d])
                monthNames.append(monthNames_d[d])
        else:
            msgwd = ''
            monthNames = []
            reserve = []
            for x in range(4, 0 ,-1):
                monthNames_x = glob.glob('%s/MSG%d/%d/%02d/*/*.nc' %(od, x, ys, ms))
                if len(monthNames_x) != 0:
                    msgwd = msgwd + '%s' %x
                    if len(monthNames) == 0:
                        monthNames = monthNames_x
                    else:
                        for fx in monthNames_x:
                            same, same_p = compareOCAListname(fx, monthNames)  # @UnusedVariable
                            if same == 1:
                                reserve.append(fx)
                            else:
                                monthNames.append(fx)
                    
        monthNames.sort()
        n_files = len(monthNames)
        if n_files == 0:
            print('No OCA Files')
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

        msgwd = int(msgwd[::-1])
#         resDict = {'lat': '', 'lon': '', \
#                    'ctp_m': '', 'ctp_v': '', 'ctp_n': '', \
#                    'cprob_m': '', 'cprob_v': '', 'cprob_n': ''}
#         for t in range(0,24):
#             resDict.update({})
        i = -1
        for ocaname in monthNames:
            print('Oca file %d of %d' %(i+1, n_files))
#             if not  (os.path.basename(ocaname) == 'W_XX-EUMETSAT-Darmstadt,OCA+,MET09+SEVIRI_C_EUMG_20070414110000_1_OR_FES_E0000_0100.nc'):
#                 continue
            i = i + 1
            if i == 0:
                tot_time = toc = tic = 0
            else:
                tot_time = tot_time + (toc - tic)
                print('Estimated time left = %d min' %(int(((tot_time / i) * (n_files - i)) /60.)))
            tic = time.time()
            #: Read OCA File
            lat, lon, sc, ctp, cprob = readOCA(ocaname)
            #: Corrupt OCA File
            if isinstance(lat, int) and (lat in [-1, -2]):
                corruptFiles.append(ocaname)
                #: Could exist as a reserve file
                same, same_p = compareOCAListname(ocaname, reserve)
                if same == 1:
                    lat, lon, sc, ctp, cprob = readOCA(reserve[same_p])
                    if isinstance(lat, int) and (lat in [-1, -2]):
                        corruptFiles.append(reserve[same_p])
                        toc = tic
                        continue
                else:
                    toc = tic
                    continue
            #: First read OCA file
            #: Nothing to add to
            if i == 0:
                ctp_v = ctp.copy()
                ctp_n = np.zeros(ctp.shape)
                ctp_n[~np.isnan(ctp)] = 1
                
                cprob_v = cprob.copy()
                cprob_n = np.zeros(cprob.shape)
                cprob_n[~np.isnan(cprob)] = 1
                
                sc_0 = (sc == 0).astype('int') #: 0=Clear
                sc_1 = (sc == 1).astype('int') #: 1=SL liquid water
                sc_2 = (sc == 2).astype('int') #: 2=SL ice water
                sc_3 = (sc == 3).astype('int') #: 3=2L
                sc_r = (~((sc == 0) | (sc == 1) | (sc == 2) | (sc == 3))).astype('int') #: nan + 10=failed retrieval
            #: Add to previus
            else:
                ctp_v = np.nansum([ctp_v, ctp], axis=0)
                ctp_n[~np.isnan(ctp)] = ctp_n[~np.isnan(ctp)] + 1
                
                cprob_v = np.nansum([cprob_v, cprob], axis=0)
                cprob_n[~np.isnan(cprob)] = cprob_n[~np.isnan(cprob)] + 1

                sc_0 = sc_0  + (sc == 0).astype('int') #: 0=Clear
                sc_1 = sc_1  + (sc == 1).astype('int') #: 1=SL liquid water
                sc_2 = sc_2  + (sc == 2).astype('int') #: 2=SL ice water
                sc_3 = sc_3  + (sc == 3).astype('int') #: 3=2L
                sc_r = sc_r  + (~((sc == 0) | (sc == 1) | (sc == 2) | (sc == 3))).astype('int') #: nan + 10=failed retrieval
            toc = time.time()

#             if i == 5:
#                 break
        
        #: Calculate Mean
        ctp_n[ctp_n==0] = np.nan
        ctp_m = ctp_v / ctp_n
        cprob_n[cprob_n==0] = np.nan
        cprob_m = cprob_v / cprob_n
        
        cfc_clo = sc_1 + sc_2 + sc_3
        cfc_tot = (cfc_clo + sc_0) * 1. #: Turn into float
        cfc_tot[cfc_tot == 0] = np.nan
        cfc = cfc_clo / cfc_tot
        resDict = {'lat': lat, 'lon': lon, \
                   'ctp_m': ctp_m, 'ctp_v': ctp_v, 'ctp_n': ctp_n, \
                   'cprob_m': cprob_m, 'cprob_v': cprob_v, 'cprob_n': cprob_n, \
                   'sc_0': sc_0, 'sc_1': sc_1, 'sc_2': sc_2, 'sc_3': sc_3, 'sc_r': sc_r, 'cfc': cfc}
        
        
        print('Total OCA read time = %d min' %(int((tot_time + (toc - tic)) / 60.)))
        print('')
        print('Save Temp file')
        tics = time.time()
        h5f = h5py.File(tempname, 'w')
        for arname, value in resDict.items():
            h5f.create_dataset('%s' %arname, data = value)
        h5f.create_dataset('MSG Name', data = [msgwd])
        h5f.create_dataset('Number of Files', data = [n_files])
        h5f.create_dataset('Corrupt Files', data = corruptFiles)
        h5f.create_dataset('Extra Files', data = reserve)
        
        h5f.close()
        tocs = time.time()
        print('Total Temp save time = %d s' %(int(tocs - tics)))
        print('')
#         if cm == True:
#             sys.exit()
        
    else:
        #: Read Temp file
        h5f = h5py.File(tempname, 'r')
        lat = h5f['lat'][:]
        lon = h5f['lon'][:]
        ctp_m = h5f['ctp_m'][:]
        cprob_m = h5f['cprob_m'][:]
        cfc = h5f['cfc'][:]
        msgwd = h5f['MSG Name'][:][0]
        n_files = h5f['Number of Files'][:][0]
        corruptFiles = h5f['Corrupt Files'][:]
        ef = h5f['Extra Files'][:]  # @UnusedVariable
        h5f.close()
    
    extra = {'MSG Name': msgwd, 'Number of Files': n_files, 'Corrupt Files': corruptFiles}
    
    return lat, lon, ctp_m, cprob_m, cfc, extra


def readOCA(fname):
    #: Some OCA files are corrupt
    try:
        ncf = nc.Dataset(fname)
    except OSError:
        return -1, -1, -1, -1, -1
    except:
        return -2, -2, -2, -2, -2
    
    ncf.variables.keys()
    lat = ncf.variables['latitude'][:].data.astype(np.float64)
    lat[ncf.variables['latitude'][:].mask] = np.nan
    
    lon = ncf.variables['longitude'][:].data.astype(np.float64)
    lon[ncf.variables['longitude'][:].mask] = np.nan

    sc = ncf.variables['scene_classification'][:].data.astype(np.float64)
    sc[ncf.variables['scene_classification'][:].mask] = np.nan
    ctp = ncf.variables['cloud_top_pressure'][:].data.astype(np.float64)
    ctp[ncf.variables['cloud_top_pressure'][:].mask] = np.nan
#     ctpll = ncf.variables['cloud_top_pressure_lower_layer'][:].data.astype(np.float64)
#     ctpll[ncf.variables['cloud_top_pressure_lower_layer'][:].mask] = np.nan
    cprob = ncf.variables['cloud_probability'][:].data.astype(np.float64)
    cprob[ncf.variables['cloud_probability'][:].mask] = np.nan
    ncf.close()
    return lat, lon, sc, ctp, cprob


def readClas(fname, typ='cfc'):
    ncf = nc.Dataset(fname)
    lat = ncf.variables['lat'][:].data.astype(np.float64)
    if ncf.variables['lat'][:].mask.ndim != 0:
        lat[ncf.variables['lat'][:].mask] = np.nan
    lat_b = ncf.variables['lat_bnds'][:].data.astype(np.float64)
    if ncf.variables['lat_bnds'][:].mask.ndim != 0:
        lat_b[ncf.variables['lat_bnds'][:].mask] = np.nan
    
    lon = ncf.variables['lon'][:].data.astype(np.float64)
    if ncf.variables['lon'][:].mask.ndim != 0:
        lon[ncf.variables['lon'][:].mask] = np.nan
    lon_b = ncf.variables['lon_bnds'][:].data.astype(np.float64)
    if ncf.variables['lon_bnds'][:].mask.ndim != 0:
        lon_b[ncf.variables['lon_bnds'][:].mask] = np.nan
    
    time = ncf.variables['time'][:].data.astype(np.float64)
    if ncf.variables['time'][:].mask.ndim != 0:
        time[ncf.variables['time'][:].mask] = np.nan
    time_b = ncf.variables['time_bnds'][:].data.astype(np.float64)
    if ncf.variables['time_bnds'][:].mask.ndim != 0:
        time_b[ncf.variables['time_bnds'][:].mask] = np.nan
    
    vara = []
    if typ == 'cfc':
        vara.append(ncf.variables['cfc'])
        vara.append(ncf.variables['cma_prob'])
    elif typ == 'cto':
        vara.append(ncf.variables['ctp'])
    var = []
    for v in vara:
        dat = v[:].data.astype(np.float64)
        if v[:].mask.ndim != 0:
            dat[v[:].mask] = np.nan
        var.append(dat)
    ncf.close()
    return lat, lon, time, var, lat_b, lon_b, time_b


def readGewax(fname):
    if os.path.isfile(fname):
        h5f = h5py.File(fname, 'r')
        latm = h5f['/']['Latitude_Midpoint'][:]
        lonm = h5f['/']['Longitude_Midpoint'][:]
        ctp = h5f['/']['Cloud_Top_Pressure_TopLayer']['Cloud_Top_Pressure_Mean_TopLayer'][:]
        cft = h5f['/']['Cloud_Amount_TopLayer']['Cloud_Amount_Mean_TopLayer'][:]
        ctp = np.where(ctp==-9999, np.nan, ctp)
        cft = np.where(cft==-9999, np.nan, cft)
        h5f.close()
    else:
        latm = lonm = ctp = cft = np.nan
    return latm, lonm, ctp, cft



def resampleProjection(dlon, dlat, slon, slat, swath_data):
    dest_def =  SwathDefinition(dlon, dlat)
    swath_def = SwathDefinition(slon, slat)
    result = resample_nearest(swath_def, swath_data, dest_def,
                              radius_of_influence=5000, fill_value=None)
    return result



if __name__ == '__main__':
    
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-l","--lat-bound", type=int, default=91,  
                        help="Latitude Boundary. Default=90")
#     parser.add_argument("-c1","--clim-max", type=int, default=500, 
#                         help="Max clim. Default=500")
#     parser.add_argument("-fn","--method-name", type=str, default='All', 
#                         help="name of the method used. Default=All")
#     parser.add_argument("-fm","--method-flag", action="store_true", default=False, 
#                         help="Do not show figure. Default=True (i.e. show figure)")
    parser.add_argument("-s","--show", action="store_false", default=True, 
                        help="Do not show figure. Default=True (i.e. show figure)")
    
    parser.add_argument("-u","--useOpt", type=int, default=10,  
                        help="use option. Default=10")
    parser.add_argument("-y","--year", type=int, default=0,  
                        help="Year. Default=2004-2019")
    parser.add_argument("-m","--month", type=int, default=0,  
                        help="Month. Default=0")
    parser.add_argument("-o","--om", action="store_true", default=False, 
                        help="Create monthly mean of OCA")
    
    
    args = parser.parse_args()
    show = args.show
    use_opt = args.useOpt
    
    
    clasMainDir = '/home/foua/data_links/data/cloud_products/CLAAS-3/L3'
    ocaMainDir_msg = '/home/foua/data_links/data/cloud_products/OCA/Release1'
    
    
    
    gewexMainDir = '/nobackup/smhid19/proj/foua/data/satellite/calipso_monthly_mean_GEWEX/h5_converted'
    plotDir = 'Plots/Compare_OCA/Test'
    tempfiles = '/TempFiles/OCAMean'
    if args.year == 0:
        years = range(2004,2021)
    else:
        years = [args.year]
    if args.month == 0:
        months = range(1,13)
    else:
        months = [args.month]
    day = 1
    
    # cfc = cloud fraction
    # cph = cloud types  
    # cto = cloud top hPa
    # cwp  jch
    
    lonremap = np.asarray([*range(-180, 180)]) + 0.5
    latremap = np.asarray([*range(-90, 90)]) + 0.5
    lonremap_mesh, latremap_mesh = np.meshgrid(lonremap, latremap)
    
    yearmonths = []
    msgNums = []
    plotMaps_done = False
    clas_ctp = np.zeros(len(months) * len(years))
    oca_ctp = np.zeros(len(months) * len(years))
    gewex_ctp = np.zeros(len(months) * len(years))
    clas_cfc = np.zeros(len(months) * len(years))
    oca_cfc = np.zeros(len(months) * len(years))
    gewex_cfc = np.zeros(len(months) * len(years))
    clas_cprob50 = np.zeros(len(months) * len(years))
    oca_cprob50 = np.zeros(len(months) * len(years))
    oca_satnum = np.zeros(len(months) * len(years))
    plot_x_labels = []
    plot_x_ticks = []
    i = -1
    for year in years:
        for mon in months:
            i = i + 1
            print('%d-%02d' %(year, mon))
            yearmonths.append('%d-%02d' %(year, mon))
            
            if mon in [6]:
                plot_x_labels.append('%d-%02d' %(year, mon))
                #: TODO: Check if plot inds start at 0 or 1
                plot_x_ticks.append(i + 1)
            if (year == 2010) and (mon == 12):
                plotMaps = True
            else:
                plotMaps = False
            cfcclasDir = '%s/cfc/%d/%02d' %(clasMainDir, year, mon)
            ctoclasDir = '%s/cto/%d/%02d' %(clasMainDir, year, mon)
    #         ocaDir = '%s/%d/%02d/%02d' %(ocaMainDir, year, mon, day)
            gewexDir = '%s/%d' %(gewexMainDir, year)
        

            #: OCA
            # msgn = getMSGNum(year, mon)
            # if msgn == 0:
            #     print('No OCA data for this month')
            #     sys.exit()
            #     ocaMainDir = '%s/%s' %(ocaMainDir_msg, msg)
            #     ocaDir = '%s/%d/%02d' %(ocaMainDir, year, mon)
            olat, olon, octp, ocprob_m, ocfc, oextra = getMeanOca(ocaMainDir_msg, year, mon, args.om) #: Pa
            
            lonremap = np.asarray([*range(-180, 180)]) + 0.5
            latremap = np.asarray([*range(-90, 90)]) + 0.5
            lonremap_mesh, latremap_mesh = np.meshgrid(lonremap, latremap)
            
            if use_opt in [10, 11, 12, 13]:
                #: hard coded from OCA 2019-08
                maxLon = 81.20243072509766
                minLon = -81.20243072509766
                maxLat = 81.26285552978516
                minLat = -81.26285552978516
            elif use_opt in [20, 21, 22, 23]:
                maxLon = 50
                minLon = -50
                maxLat = 50
                minLat = -50
            
            if (not isinstance(olat, np.ndarray)) and np.isnan(olat):
                oca_ctp[i] = np.nan
                oca_cfc[i] = np.nan
#                 oca_cprob50[i] = np.nan
                msgNums.append(np.nan)
                #: hard coded from OCA 2019-08
            else:
                #: OCA Lon starts from high to the left and get smaller
                #: Other data do opposite 
                #: Hence flip OCA
                olat = olat[:, ::-1]
                olon = olon[:, ::-1]
                octp = octp[:, ::-1]
                ocfc = ocfc[:, ::-1]
#                 ocprob_m = ocprob_m[:, ::-1]
                if use_opt in [10, 11, 20, 21]:
                    ocfc = ocfc * 100 #: Procent
                    octp = octp / 100. #: Change to hPa
                elif use_opt in [12, 13, 22, 23]:
                    ocfc_masked = resampleProjection(lonremap_mesh, latremap_mesh, olon, olat, ocfc * 100)
                    octp_masked = resampleProjection(lonremap_mesh, latremap_mesh, olon, olat, octp / 100.)
                    ocfc =  np.where(ocfc_masked.mask, np.nan, ocfc_masked.data)
                    octp =  np.where(octp_masked.mask, np.nan, octp_masked.data)
                    
                    olon = lonremap_mesh
                    olat = latremap_mesh
                msgNums.append(oextra['MSG Name'])
                #: OCA
#                 maxLon = np.nanmax(olon)
#                 minLon = np.nanmin(olon)
#                 maxLat = np.nanmax(olat)
#                 minLat = np.nanmin(olat)
                ocaInd = (olon <= maxLon) & (olon >= minLon) & (olat <= maxLat) & (olat >= minLat)
#                 oca_nanind = ~np.isnan(ocprob_m)
                if use_opt in [10, 20, 12, 22]:
                    oca_ctp[i] = np.nanmean(octp[ocaInd])
                    oca_cfc[i] = np.nanmean(ocfc[ocaInd])
                elif use_opt in [11, 21, 13, 23]:
                    oca_ctp[i] = np.nanmean(np.where(~np.isnan(octp), octp * np.cos(np.deg2rad(olat)), np.nan)[ocaInd])
                    oca_cfc[i] = np.nanmean(np.where(~np.isnan(ocfc), ocfc * np.cos(np.deg2rad(olat)), np.nan)[ocaInd])
#                 oca_cprob50[i] = (ocprob_m[oca_nanind]>50).sum() / oca_nanind.sum()

            #: Clas
            cfcClasname = glob.glob('%s/CFCmm*.nc' %cfcclasDir)[0] #: md = Diurnal cycle, mm = monthly mean
            ctoClasname = glob.glob('%s/CTOmm*.nc' %ctoclasDir)[0] #: md = Diurnal cycle, mm = monthly mean, mh = ?
            cfclat, cfclon, cfctime, cfc_var, cfclat_b, cfclon_b, cfctime_b = readClas(cfcClasname, 'cfc')
            ctolat, ctolon, ctotime, ctp_var, ctolat_b, ctolon_b, ctotime_b = readClas(ctoClasname, 'cto') #: hPa
            if np.all(ctp_var[0] == 0) or ((not isinstance(cfclat, np.ndarray)) and np.isnan(cfclat)):
                clas_ctp[i] = np.nan
                clas_cfc[i] = np.nan
#                 clas_cprob50[i] = np.nan
            else:
                ctolon_mesh, ctolat_mesh = np.meshgrid(ctolon, ctolat)
                if use_opt in [10, 11, 20, 21]:
                    c_ctp = ctp_var[0][0]
                    c_cfc = cfc_var[0][0]
                elif use_opt in [12, 13, 22, 23]:
#                     c_ctp = resampleProjection(lonremap_mesh, latremap_mesh, ctolon_mesh, ctolat_mesh, ctp_var[0][0])
#                     c_cfc = resampleProjection(lonremap_mesh, latremap_mesh, ctolon_mesh, ctolat_mesh, cfc_var[0][0])
#                     ctolon = lonremap
#                     ctolat = latremap
#                     ctolon_mesh = lonremap_mesh
#                     ctolat_mesh = latremap_mesh
                    c_ctp = ctp_var[0][0][10::20, 10::20]
                    c_cfc = cfc_var[0][0][10::20, 10::20]
                    ctolon = ctolon[10::20]
                    ctolat = ctolat[10::20]
                    ctolon_mesh = ctolon_mesh[10::20, 10::20]
                    ctolat_mesh = ctolat_mesh[10::20, 10::20]
#                 c_cprob = cfc_var[1]
                #: TODO: Right order?Update. I think it is
                #: It wasent I think
                #: Clas
                clasInd = (ctolon_mesh <= maxLon) & (ctolon_mesh >= minLon) & (ctolat_mesh <= maxLat) & (ctolat_mesh >= minLat)
                if use_opt in [10, 20, 12, 22]:
                    clas_ctp[i] = np.nanmean(c_ctp[clasInd])
                    clas_cfc[i] = np.nanmean(c_cfc[clasInd])
                elif use_opt in [11, 21, 13, 23]:
                    clas_ctp[i] = np.nanmean(np.where(~np.isnan(c_ctp), c_ctp *  np.cos(np.deg2rad(ctolat_mesh)), np.nan)[clasInd])
                    clas_cfc[i] = np.nanmean(np.where(~np.isnan(c_cfc), c_cfc *  np.cos(np.deg2rad(ctolat_mesh)), np.nan)[clasInd])
#                 clas_cprob50[i] = (c_cprob[clas_nanind]>50).sum() / clas_nanind.sum()
            
            #: Gewax
            gewexname = '%s/CAL_LID_L3_GEWEX_Cloud-Standard-V1-00.%i-%02dA.h5' %(gewexDir, year, mon)
            glat, glon, gctp, gcft = readGewax(gewexname) #: hPa
            if (not isinstance(glat, np.ndarray)) and np.isnan(glat):
                gewex_ctp[i] = np.nan
                gewex_cfc[i] = np.nan
            else:
                glon_mesh, glat_mesh = np.meshgrid(glon, glat)
                #: gewex
                gewexInd = (glon_mesh <= maxLon) & (glon_mesh >= minLon) & (glat_mesh <= maxLat) & (glat_mesh >= minLat)
                gcft = gcft * 100.
                if use_opt in [10, 20, 12, 22]:
                    gewex_ctp[i] = np.nanmean(gctp[gewexInd])
                    gewex_cfc[i] = np.nanmean(gcft[gewexInd])
                elif use_opt in [11, 21, 13, 23]:
                    gewex_ctp[i] = np.nanmean(np.where(~np.isnan(gctp),  gctp *  np.cos(np.deg2rad(glat_mesh)), np.nan)[gewexInd])
                    gewex_cfc[i] = np.nanmean(np.where(~np.isnan(gcft),  gcft *  np.cos(np.deg2rad(glat_mesh)), np.nan)[gewexInd])
            
            
            
            
            
            if plotMaps:
                ocaMap_ctp = np.where(ocaInd, octp, np.nan)
                ocaMap_cfc = np.where(ocaInd, ocfc, np.nan)
                ocaExtent = (np.nanmax(olon), np.nanmin(olon), np.nanmax(olat), np.nanmin(olat))
                claasMap_ctp = np.where(clasInd, c_ctp, np.nan)
                claasMap_cfc = np.where(clasInd, c_cfc, np.nan)
                claasExtent = (np.nanmax(ctolon), np.nanmin(ctolon), np.nanmax(ctolat), np.nanmin(ctolat))
                gewexMap_ctp = np.where(gewexInd, gctp, np.nan)
                gewexMap_cfc = np.where(gewexInd, gcft, np.nan)
                gewexExtent = (np.nanmax(glon), np.nanmin(glon), np.nanmax(glat), np.nanmin(glat))
                plotMaps_done = True
#     class cartopy.crs.Orthographic(central_longitude=0.0, central_latitude=0.0, globe=None)[source]
# 
# 
# 
    
#     octp_use = []
#     cctp_use = []
#     for i in range(olon.shape[0]-1):
    
    oca_ctp_plot = oca_ctp
    clas_ctp_plot = clas_ctp
    gewex_ctp_plot = gewex_ctp
    oca_cfc_plot = oca_cfc
    clas_cfc_plot = clas_cfc
    gewex_cfc_plot = gewex_cfc
    msg_colour = {1: 'g', 2: 'olivedrab', 3: 'lime', 4: 'darkgreen', 12: 'r', 13: 'violet', 23: 'hotpink', 34: 'darkviolet'}
    plot_x = [*range(1, len(clas_ctp_plot) + 1)]
    if args.year == 0:
        figend = ''
    else:
        figend = '-%s' %year
    figend = figend + '_opt-%d' %use_opt
    title_extra = '\n Area = +-%d' %(int(maxLon))

    if use_opt in [10, 20, 12, 22]:
        title_extra = title_extra + ', cos(lat) = False'
    elif use_opt in [11, 21, 13, 23]:
        title_extra = title_extra + ', cos(lat) = True'
    if use_opt in [10, 11, 20, 21]:
        title_extra = title_extra + ', Resample = False'
    elif use_opt in [12, 13, 22, 23]:
        title_extra = title_extra + ', Resample = True'
        
    #: Only plot time series if there are enough months
    if len(months) >= 12:
        #: To make sure label for OCA is only used once
        msg_label_used = {1: False, 2: False, 3: False, 4: False, 12: False, 13: False, 23: False, 34: False}
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(1,1,1)
        ax.plot(plot_x, clas_ctp_plot, label='CLAAS')
        ax.plot(plot_x, gewex_ctp_plot, label='GEWEX')
#         ax.plot(plot_x, oca_ctp, label='OCA', color='k')
        #: Plot different MSG (OCA)
        for i in range(len(clas_ctp_plot)-1):
            msgN = msgNums[i]
            if np.isnan(msgN):
                continue
            msgN_string = '%d' %msgN
            if len(msgN_string) == 1:
                msgLabel = 'MSG%s' %msgN_string
            elif len(msgN_string) == 2:
                msgLabel = 'MSG%s/%s' %(msgN_string[0], msgN_string[1])
            else:
                print('hmm MSG satellites in Plot')
                pdb.set_trace()
            if not msg_label_used[msgN]:
                ax.plot([plot_x[i], plot_x[i+1]], [oca_ctp_plot[i], oca_ctp_plot[i + 1]], color=msg_colour[msgN], linestyle='-',  label='OCA - %s' %msgLabel)
                msg_label_used[msgN] = True
            else:
                ax.plot([plot_x[i], plot_x[i+1]], [oca_ctp_plot[i], oca_ctp_plot[i + 1]], color=msg_colour[msgN], linestyle='-')
        
        
    #===========================================================================
    # for s in set(msgNums):
    #     if np.isnan(s):
    #         continue
    #     else:
    #         ocapltind = (msgNums == s)
    #         plot_y_oca = np.where(msgNums == s, oca_ctp, np.nan)
    #         s_string = '%d' %s
    #         if len(s_string) == 1:
    #             msgLabel = 'MSG%s' %s_string
    #         elif len(s_string) == 2:
    #             msgLabel = 'MSG%s/%s' %(s_string[0], s_string[1])
    #         else:
    #             print('hmm MSG satellites in Plot')
    #             pdb.set_trace()
    #     plot_y_oca = np.where(msgNums == s, oca_ctp, np.nan)
    #===========================================================================
            
        ax.set_xlabel('Month')
    #     ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    #     ax.set_xticklabels(['2010-01', '2010-02', '2010-03', '2010-04', '2010-05', '2010-06', '2010-07', '2010-08', '2010-09', '2010-10', '2010-11', '2010-12'])
        ax.set_xticks(plot_x_ticks)
        ax.set_xticklabels(plot_x_labels)
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_ylabel('Cloud Top Pressure [hPa]')
        ax.legend(fontsize=8, bbox_to_anchor=(1.02, 1))
        ax.set_title('Monthly Mean Cloud Top Pressure' + title_extra)
        fig.autofmt_xdate()
        plt.tight_layout()
        figname = '%s/ctp_month-mean%s' %(plotDir, figend)
        fig.savefig(figname +'.png')
        print(figname)
    
        msg_label_used = {1: False, 2: False, 3: False, 4: False, 12: False, 13: False, 23: False, 34: False}
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(1,1,1)
    #     ax.plot(months, oca_cfc, label='OCA')
        ax.plot(plot_x, clas_cfc_plot, label='CLAAS')
        ax.plot(plot_x, gewex_cfc_plot, label='GEWEX')
        for i in range(len(clas_cfc_plot)-1):
            msgN = msgNums[i]
            if np.isnan(msgN):
                continue
            msgN_string = '%d' %msgN
            if len(msgN_string) == 1:
                msgLabel = 'MSG%s' %msgN_string
            elif len(msgN_string) == 2:
                msgLabel = 'MSG%s/%s' %(msgN_string[0], msgN_string[1])
            else:
                print('hmm MSG satellites in Plot')
                pdb.set_trace()
            if not msg_label_used[msgN]:
                ax.plot([plot_x[i], plot_x[i+1]], [oca_cfc_plot[i], oca_cfc_plot[i + 1]], color=msg_colour[msgN], linestyle='-',  label='OCA - %s' %msgLabel)
                msg_label_used[msgN] = True
            else:
                ax.plot([plot_x[i], plot_x[i+1]], [oca_cfc_plot[i], oca_cfc_plot[i + 1]], color=msg_colour[msgN], linestyle='-')
        ax.set_xlabel('Month')
        ax.set_xticks(plot_x_ticks)
        ax.set_xticklabels(plot_x_labels)
    #     ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_ylabel('Cloud Fraction [%]')
        ax.set_ylim([40, 75])
        ax.legend(fontsize=8, bbox_to_anchor=(1.02, 1))
    
        ax.set_title('Monthly Mean Cloud Fraction' + title_extra)
        fig.autofmt_xdate()
        plt.tight_layout()
        figname = '%s/cfc_month-mean%s' %(plotDir, figend)
        fig.savefig(figname +'.png')
        print(figname)
    
    if False:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(months, oca_cprob50, label='OCA')
        ax.plot(months, clas_cprob50, label='CLAAS')
    #     ax.plot(months, gewex_cfc, label='GEWEX')
        ax.set_xlabel('Month')
        ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
        ax.set_xticklabels(['2010-01', '2010-02', '2010-03', '2010-04', '2010-05', '2010-06', '2010-07', '2010-08', '2010-09', '2010-10', '2010-11', '2010-12'])
    #     ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_ylabel('Cloud Prob > 50% [%]')
        ax.legend(fontsize=8)
        ax.set_title('Monthly Mean Cloud Prob > 50% for 2010')
        fig.autofmt_xdate()
        plt.tight_layout()
        figname = '%s/cprob50_month-mean-%s' %(plotDir, year) + '_opt-%d' %use_opt
        fig.savefig(figname +'.png')

    #: Plot maps
    if plotMaps_done:
        oca_proj = ccrs.Geostationary(satellite_height=35786000)
        #: CTP
        plt.figure(figsize=(3, 3))
        fig = plt.figure()
        ax = fig.add_subplot(111,projection=ccrs.Orthographic())
    #     ax.gridlines()
        im = ax.imshow(ocaMap_ctp[:,::-1], extent = ocaExtent, transform=ccrs.PlateCarree(), vmin=100, vmax=1000)#, alpha=0.5)
#         im = ax.imshow(ocaMap_ctp, origin ='lower', extent = ocaExtent, transform=oca_proj, vmin=100, vmax=1000)#, alpha=0.5)
        ax.coastlines(resolution='110m')#, color='r')
        ax.set_title('OCA, %d-%02d' %(year, mon) + title_extra)
        cbar = fig.colorbar(im)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
        cbar.set_label('Cloud Top Pressure [hPa]')
        figname = '%s/map_ctp_oca_month-mean-%d-%02d' %(plotDir, year, mon) + '_opt-%d' %use_opt
        fig.savefig(figname +'.png')
        print(figname)
        
        plt.figure(figsize=(3, 3))
        fig = plt.figure()
        ax = fig.add_subplot(111,projection=ccrs.Orthographic())
    #     ax.gridlines()
        im = ax.imshow(claasMap_ctp[:,::-1], extent = claasExtent, transform=ccrs.PlateCarree(), vmin=100, vmax=1000)#, alpha=0.5)
        ax.coastlines(resolution='110m')#, color='r')
        ax.set_title('CLAAS, %d-%02d' %(year, mon) + title_extra)
        cbar = fig.colorbar(im)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
        cbar.set_label('Cloud Top Pressure [hPa]')
        figname = '%s/map_ctp_claas_month-mean-%d-%02d' %(plotDir, year, mon) + '_opt-%d' %use_opt
        fig.savefig(figname +'.png')
        print(figname)
        
        plt.figure(figsize=(3, 3))
        fig = plt.figure()
        ax = fig.add_subplot(111,projection=ccrs.Orthographic())#central_longitude=-80))
    #     ax.gridlines()
        im = ax.imshow(gewexMap_ctp[:,::-1], extent = gewexExtent, transform=ccrs.PlateCarree(), vmin=100, vmax=1000)#, alpha=0.5)
    #     ax.imshow(g_ctp_ocaLL, origin='lower', extent = (maxLon, minLon, maxLat, minLat), transform=ccrs.PlateCarree())#, alpha=0.5)
        ax.coastlines(resolution='110m')#, color='r')
        ax.set_title('GEWEX, %d-%02d' %(year, mon) + title_extra)
        cbar = fig.colorbar(im)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
        cbar.set_label('Cloud Top Pressure [hPa]')
        figname = '%s/map_ctp_gewex_month-mean-%d-%02d' %(plotDir, year, mon) + '_opt-%d' %use_opt
        fig.savefig(figname +'.png')
        print(figname)
        
        #: CFC
        
        
#         plt.figure(figsize=(3, 3))
#         fig = plt.figure()
#         ax = fig.add_subplot(111,projection=ccrs.Orthographic())
#     #     ax.gridlines()
#         im = ax.imshow(ocfc_re, origin ='lower', extent = gewexExtent, transform=ccrs.PlateCarree(), vmin=0, vmax=100)#, alpha=0.5)
#         ax.coastlines(resolution='110m')#, color='r')
#         ax.set_title('OCA, %d-%02d' %(year, mon) + title_extra)
#         cbar = fig.colorbar(im)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
#         cbar.set_label('Cloud Fraction [%]')
#         figname = 'test'
#         fig.savefig(figname +'.png')
#         print(figname)
        
        
        
        plt.figure(figsize=(3, 3))
        fig = plt.figure()
        ax = fig.add_subplot(111,projection=ccrs.Orthographic())
    #     ax.gridlines()
        im = ax.imshow(ocaMap_cfc[:,::-1], extent = ocaExtent, vmin=0, vmax=100, transform=ccrs.PlateCarree())#, alpha=0.5)
        ax.coastlines(resolution='110m')#, color='r')
        ax.set_title('OCA, %d-%02d' %(year, mon) + title_extra)
        cbar = fig.colorbar(im)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
        cbar.set_label('Cloud Fraction [%]')
        figname = '%s/map_cfc_oca_month-mean-%d-%02d' %(plotDir, year, mon) + '_opt-%d' %use_opt
        fig.savefig(figname +'.png')
        print(figname)
        
        plt.figure(figsize=(3, 3))
        fig = plt.figure()
        ax = fig.add_subplot(111,projection=ccrs.Orthographic())
    #     ax.gridlines()
        im = ax.imshow(claasMap_cfc[:,::-1], extent = claasExtent, transform=ccrs.PlateCarree(), vmin=0, vmax=100)#, alpha=0.5)
        ax.coastlines(resolution='110m')#, color='r')
        ax.set_title('CLAAS, %d-%02d' %(year, mon) + title_extra)
        cbar = fig.colorbar(im)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
        cbar.set_label('Cloud Fraction [%]')
        figname = '%s/map_cfc_claas_month-mean-%d-%02d' %(plotDir, year, mon) + '_opt-%d' %use_opt
        fig.savefig(figname +'.png')
        print(figname)
        
        plt.figure(figsize=(3, 3))
        fig = plt.figure()
        ax = fig.add_subplot(111,projection=ccrs.Orthographic())#central_longitude=-80))
    #     ax.gridlines()
        im = ax.imshow(gewexMap_cfc[:,::-1], extent = gewexExtent, transform=ccrs.PlateCarree(), vmin=0, vmax=100)#, alpha=0.5)
    #     ax.imshow(g_ctp_ocaLL, origin='lower', extent = (maxLon, minLon, maxLat, minLat), transform=ccrs.PlateCarree())#, alpha=0.5)
        ax.coastlines(resolution='110m')#, color='r')
        ax.set_title('GEWEX, %d-%02d' %(year, mon) + title_extra)
        cbar = fig.colorbar(im)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
        cbar.set_label('Cloud Fraction [%]')
        figname = '%s/map_cfc_gewex_month-mean-%d-%02d' %(plotDir, year, mon) + '_opt-%d' %use_opt
        fig.savefig(figname +'.png')    
        print(figname)





    
    pdb.set_trace()
    
    
    
    
    
    