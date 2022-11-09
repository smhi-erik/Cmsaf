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
                    

def getMeanOca(od, ys, ms, tid=-1, cm=False, lt=True):
    #: od = dir
    #: ys = year
    #: ms = month
    #: cm = create monthly mean
    #: lt = load temp files
    
    
    tempname = os.path.join('TempFiles', 'OCAMean', 'oca_monthly-mean_%d-%02d.h5' %(ys, ms))
    if tid >= 0:
        tempname = tempname.replace('.h5', '_t-%02d.h5' %tid)
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
    
    #: Controle for corrupt files. 
    #: OSError is usually that it can not be open
    #: KeyError data is missing
    #: Only d this when creating new month files i.e. cm=True
    important_variables = ['lat', 'lon', 'ctp_m', 'cprob_m', 'cfc', 'MSG Name', 'Number of Files', 'Corrupt Files', 'Extra Files']
    if cm and os.path.isfile(tempname):
        try:
            h5f = h5py.File(tempname, 'r')
            for iv in important_variables:
                temp = h5f[iv]  # @UnusedVariable
        except (OSError, KeyError):
            os.remove(tempname)
        else:
            h5f.close()

    if (lt == False) or (not os.path.isfile(tempname)):
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
        #: Remove files not used in time
        monthNames_tid = []
        if tid >= 0:
            #===================================================================
            # ncf = nc.Dataset('/home/foua/data_links/data/cloud_products/CLAAS-3/L3/cfc/2010/01/CFCmd20100101000000419SVMSG01MA.nc')
            # claas_time_bnds = (ncf.variables['time_bnds'][:] - ncf.variables['time_bnds'][:][0,0]) * 24
            # ncf.close()
            #===================================================================
            #: Def from claas, see above
            tb_min = tid
            tb_max = tid + 1
            for fn in monthNames:
                oca_hour = int(os.path.basename(fn).split('_')[4][8:10])
                if (oca_hour >= tb_min) and (oca_hour < tb_max):
                    monthNames_tid.append(fn)
            monthNames = monthNames_tid
        
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
        print(tempname)
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


def addExtent(obt, typ, lon1d, lon2d):
#: Add lon extent to make sure it gets the same size as the other two
    if obt.ndim == 1:
        retv = np.zeros(lon1d.shape)
    else:
        retv = np.zeros(lon2d.shape)
    
    if typ == 'lon':
        if obt.ndim == 1:
            #: Add before
            #: Add -180 - -90  by adding 90 on first part
            retv[0:90] = obt[0:90] - 90
            #: Add middle
            #: Add -90 - 90
            #: Same as org data
            retv[90:270] = obt
            #: Add After
            #: Add 90 - 180 by adding 90 onlast part
            retv[270:360] = obt[90:] + 90
        else:
            #: Same but for 2D
            retv[:, 0:90] = obt[:, 0:90] - 90
            retv[:, 90:270] = obt
            retv[:, 270:360] = obt[:, 90:] + 90
    elif typ == 'lat':
        #: Lat is same in this axes so just add the same row to all
        retv[:, 0:90] = obt[:, 0:90]
        retv[:, 90:270] = obt
        retv[:, 270:360] = obt[:, 90:]
    elif typ == 'data':
        #: No data except org
        #: Everything before and anfter org should be nan
        retv = retv + np.nan
        #: Middle is same as org
        retv[:, 90:270] = obt
    return retv

if __name__ == '__main__':
    
    
    import argparse
    parser = argparse.ArgumentParser()
#     parser.add_argument("-l","--lat-bound", type=int, default=91,  
#                         help="Latitude Boundary. Default=90")
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
    parser.add_argument("-l","--lt", action="store_false", default=True, 
                        help="Create monthly mean of OCA")
    parser.add_argument("-t","--tid", type=int, default=-1,  
                        help="Hour of day. Default=-1 i.e. the whole day")
    
    
    args = parser.parse_args()
    show = args.show
    use_opt = args.useOpt
    
    
    clasMainDir = '/home/foua/data_links/data/cloud_products/CLAAS-3/L3'
    ocaMainDir_msg = '/home/foua/data_links/data/cloud_products/OCA/Release1'
    
    
    
    gewexMainDir = '/nobackup/smhid19/proj/foua/data/satellite/calipso_monthly_mean_GEWEX/h5_converted'
    plotDir = 'Plots/Compare_OCA'
    tempfiles = '/TempFiles/OCAMean'
    if args.year == 0:
        years = range(2004,2020)
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
            #: No OCA data after 2019-08
            if ((year == 2019) and (mon > 8)) or (year > 2019):
                break
            i = i + 1
            if args.tid == -1:
                print('%d-%02d' %(year, mon))
            else:
                print('%d-%02d-t%02d' %(year, mon, args.tid))
            yearmonths.append('%d-%02d' %(year, mon))
            
            #: Labels used for plotting time series
            if (len(years) == 1) or ((len(years) >= 1) and (mon in [6])):
                plot_x_labels.append('%d-%02d' %(year, mon))
                plot_x_ticks.append(i + 1)
            #: Year / mont that is plotted on a map
            if (year == 2010) and (mon == 12):
#             if (year == 2005) and (mon == 3):
                plotMaps = True
            else:
                plotMaps = False
            #: Dirs
            cfcclasDir = '%s/cfc/%d/%02d' %(clasMainDir, year, mon)
            ctoclasDir = '%s/cto/%d/%02d' %(clasMainDir, year, mon)
            gewexDir = '%s/%d' %(gewexMainDir, year)
        
            
            #: Read data
            #: OCA
            olat, olon, octp, ocprob_m, ocfc, oextra = getMeanOca(ocaMainDir_msg, year, mon, args.tid, args.om, args.lt) #: Pa
            if args.om:
                continue
            #: Controle if the data existed for the particular month
            if (not isinstance(olat, np.ndarray)) and np.isnan(olat):
                no_oca_data = True
            else:
                no_oca_data = False
            #: Clas
            cfcClasname = glob.glob('%s/CFCmm*.nc' %cfcclasDir)[0] #: md = Diurnal cycle, mm = monthly mean
            ctoClasname = glob.glob('%s/CTOmm*.nc' %ctoclasDir)[0] #: md = Diurnal cycle, mm = monthly mean, mh = ?
            cfclat, cfclon, cfctime, cfc_var, cfclat_b, cfclon_b, cfctime_b = readClas(cfcClasname, 'cfc')
            ctolat, ctolon, ctotime, ctp_var, ctolat_b, ctolon_b, ctotime_b = readClas(ctoClasname, 'cto') #: hPa
            #: Controle if the data existed for the particular month
            if np.all(ctp_var[0] == 0) or ((not isinstance(cfclat, np.ndarray)) and np.isnan(cfclat)):
                no_clas_data = True
            else:
                no_clas_data = False
            #: Gewax
            gewexname = '%s/CAL_LID_L3_GEWEX_Cloud-Standard-V1-00.%i-%02dA.h5' %(gewexDir, year, mon)
            glat, glon, gctp, gcft = readGewax(gewexname) #: hPa
            #: Controle if the data existed for the particular month
            if (not isinstance(glat, np.ndarray)) and np.isnan(glat):
                no_gewex_data = True
            else:
                no_gewex_data = False
            #: OCA
            #: Clas
            #: Gewax
            
            #: Remapp OCA and CLAS
            #: OCA
            if no_oca_data:
#                 oca_cprob50[i] = np.nan
                msgNums.append(np.nan)
                #: Create a map where there is OCA Data
                #: Only used when remapp
                #: This is a dummy, all true to not create conflicts
                if use_opt in [12, 13, 22, 23]:
                    oca_true_cfc = np.ones(lonremap_mesh.shape).astype('bool')
                    oca_true_ctp = np.ones(lonremap_mesh.shape).astype('bool')
                    oca_ind_latlon = np.ones(lonremap_mesh.shape).astype('bool')
            else:
                #: OCA Lon starts from high to the left and get smaller
                #: Other data do opposite 
                #: Hence flip OCA
                olat = olat[:, ::-1]
                olon = olon[:, ::-1]
                octp = octp[:, ::-1]
                ocfc = ocfc[:, ::-1]
                ocfc = ocfc * 100 #: Procent
                octp = octp / 100. #: Change to hPa
                if use_opt in [12, 13, 22, 23]:
                    #: Remapp
                    ocfc_masked = resampleProjection(lonremap_mesh, latremap_mesh, olon, olat, ocfc)
                    octp_masked = resampleProjection(lonremap_mesh, latremap_mesh, olon, olat, octp)
                    #: Remove mask
                    ocfc =  np.where(ocfc_masked.mask, np.nan, ocfc_masked.data)
                    octp =  np.where(octp_masked.mask, np.nan, octp_masked.data)
                    #: Use remapp lon/lat
                    olon = lonremap_mesh
                    olat = latremap_mesh
                #: What satellite (MSG) the OCA is from
                msgNums.append(oextra['MSG Name'])
                #: OCA Data
                oca_true_cfc = ~np.isnan(ocfc)
                oca_true_ctp = ~np.isnan(octp)
                #: OCA Lat/Lon
                oca_ind_latlon = (olon <= maxLon) & (olon >= minLon) & (olat <= maxLat) & (olat >= minLat)
                
            #: Clas
            if no_clas_data:
                if use_opt in [12, 13, 22, 23]:
                    clas_true_cfc = np.ones(lonremap_mesh.shape).astype('bool')
                    clas_true_ctp = np.ones(lonremap_mesh.shape).astype('bool')
                    clas_ind_latlon = np.ones(lonremap_mesh.shape).astype('bool')
            else:
                #: TODO: Right order?Update. I think it is
                clon_mesh, clat_mesh = np.meshgrid(ctolon, ctolat)
                cctp = ctp_var[0][0]
                ccfc = cfc_var[0][0]
                if use_opt in [12, 13, 22, 23]:
                    #: No remapp just pick every 20th, start 10 in
                    cctp = cctp[10::20, 10::20]
                    ccfc = ccfc[10::20, 10::20]
                    ctolon = ctolon[10::20]
                    ctolat = ctolat[10::20]
                    clon_mesh = clon_mesh[10::20, 10::20]
                    clat_mesh = clat_mesh[10::20, 10::20]
                    
                    #: Needs to be same extent
                    ctolon = addExtent(ctolon, 'lon', lonremap, lonremap_mesh)
                    clon_mesh = addExtent(clon_mesh, 'lon', lonremap, lonremap_mesh)
                    clat_mesh = addExtent(clat_mesh, 'lat', lonremap, lonremap_mesh)
                    cctp = addExtent(cctp, 'data', lonremap, lonremap_mesh)
                    ccfc = addExtent(ccfc, 'data', lonremap, lonremap_mesh)

                    
                clas_true_cfc = ~np.isnan(ccfc)
                clas_true_ctp = ~np.isnan(cctp)
                clas_ind_latlon = (clon_mesh <= maxLon) & (clon_mesh >= minLon) & (clat_mesh <= maxLat) & (clat_mesh >= minLat)
            
            #: Gewex no remapp but needs to been taken care of
            if no_gewex_data:
                if use_opt in [12, 13, 22, 23]:
                    gewex_true_cfc = np.ones(lonremap_mesh.shape).astype('bool')
                    gewex_true_ctp = np.ones(lonremap_mesh.shape).astype('bool')
                    gewex_ind_latlon = np.ones(lonremap_mesh.shape).astype('bool')
            else:
                #: Procent
                #: #: gcft is gcfc for top layer
                gcfc = gcft * 100.
                gewex_true_cfc = ~np.isnan(gcfc)
                gewex_true_ctp = ~np.isnan(gctp)
                glon_mesh, glat_mesh = np.meshgrid(glon, glat)
                gewex_ind_latlon = (glon_mesh <= maxLon) & (glon_mesh >= minLon) & (glat_mesh <= maxLat) & (glat_mesh >= minLat)
            
            
            #: Create Inds
            
            if [12, 13, 22, 23]:
                true_cfc = oca_true_cfc & clas_true_cfc & gewex_true_cfc
                true_ctp = oca_true_ctp & clas_true_ctp & gewex_true_ctp
                ind_latlon = oca_ind_latlon & clas_ind_latlon & gewex_ind_latlon
                
                ocaIndcfc = true_cfc & ind_latlon
                ocaIndctp = true_ctp & ind_latlon
                clasIndcfc = true_cfc & ind_latlon
                clasIndctp = true_ctp & ind_latlon
                gewexIndcfc = true_cfc & ind_latlon
                gewexIndctp = true_ctp & ind_latlon
            else:
                if not no_oca_data:
                    ocaIndcfc = oca_ind_latlon
                    ocaIndctp = oca_ind_latlon
                if not no_clas_data:
                    clasIndcfc = clas_ind_latlon
                    clasIndctp = true_ctp & ind_latlon
                if not no_gewex_data:
                    gewexIndcfc = true_cfc & ind_latlon
                    gewexIndctp = true_ctp & ind_latlon
                
            #: Calculate mean
            #: OCA
            if not no_oca_data:
                if use_opt in [10, 20, 12, 22]:
                    oca_ctp[i] = np.nanmean(octp[ocaIndctp])
                    oca_cfc[i] = np.nanmean(ocfc[ocaIndcfc])
                elif use_opt in [11, 21, 13, 23]:
                    oca_ctp[i] = np.nanmean(np.where(~np.isnan(octp), octp * np.cos(np.deg2rad(olat)), np.nan)[ocaIndctp])
                    oca_cfc[i] = np.nanmean(np.where(~np.isnan(ocfc), ocfc * np.cos(np.deg2rad(olat)), np.nan)[ocaIndcfc])
            else:
                oca_ctp[i] = np.nan
                oca_cfc[i] = np.nan
            #: Clas
            if not no_clas_data:
                if use_opt in [10, 20, 12, 22]:
                    clas_ctp[i] = np.nanmean(cctp[clasIndctp])
                    clas_cfc[i] = np.nanmean(ccfc[clasIndcfc])
                elif use_opt in [11, 21, 13, 23]:
                    clas_ctp[i] = np.nanmean(np.where(~np.isnan(cctp), cctp *  np.cos(np.deg2rad(clat_mesh)), np.nan)[clasIndctp])
                    clas_cfc[i] = np.nanmean(np.where(~np.isnan(ccfc), ccfc *  np.cos(np.deg2rad(clat_mesh)), np.nan)[clasIndcfc])
#                 clas_cprob50[i] = (c_cprob[clas_nanind]>50).sum() / clas_nanind.sum()
            else:
                clas_ctp[i] = np.nan
                clas_cfc[i] = np.nan
            #: Gewax
            if not no_gewex_data:
                if use_opt in [10, 20, 12, 22]:
                    gewex_ctp[i] = np.nanmean(gctp[gewexIndctp])
                    gewex_cfc[i] = np.nanmean(gcfc[gewexIndcfc])
                elif use_opt in [11, 21, 13, 23]:
                    gewex_ctp[i] = np.nanmean(np.where(~np.isnan(gctp),  gctp *  np.cos(np.deg2rad(glat_mesh)), np.nan)[gewexIndctp])
                    gewex_cfc[i] = np.nanmean(np.where(~np.isnan(gcfc),  gcfc *  np.cos(np.deg2rad(glat_mesh)), np.nan)[gewexIndcfc])
            else:
                gewex_ctp[i] = np.nan
                gewex_cfc[i] = np.nan
            if plotMaps:
#                 ocaMap_ctp = np.where(oca_ind_latlon, octp, np.nan)
#                 ocaMap_cfc = np.where(oca_ind_latlon, ocfc, np.nan)
#                 ocaExtent = (np.nanmax(olon), np.nanmin(olon), np.nanmax(olat), np.nanmin(olat))
#                 claasMap_ctp = np.where(clas_ind_latlon, cctp, np.nan)
#                 claasMap_cfc = np.where(clas_ind_latlon, ccfc, np.nan)
#                 claasExtent = (np.nanmax(ctolon), np.nanmin(ctolon), np.nanmax(ctolat), np.nanmin(ctolat))
#                 gewexMap_ctp = np.where(gewex_ind_latlon, gctp, np.nan)
#                 gewexMap_cfc = np.where(gewex_ind_latlon, gcfc, np.nan)
#                 gewexExtent = (np.nanmax(glon), np.nanmin(glon), np.nanmax(glat), np.nanmin(glat))
                if no_oca_data:
                    ocaMap_ctp = np.zeros(lonremap_mesh.shape) + np.nan
                    ocaMap_cfc = np.zeros(lonremap_mesh.shape) + np.nan
                    ocaExtent = (minLon, maxLon, minLat, maxLat)
                else:
                    ocaMap_ctp = np.where(ocaIndctp, octp, np.nan)
                    ocaMap_cfc = np.where(ocaIndcfc, ocfc, np.nan)
#                     ocaExtent = (np.nanmax(olon), np.nanmin(olon), np.nanmax(olat), np.nanmin(olat))
                    ocaExtent = (np.nanmin(olon), np.nanmax(olon), np.nanmin(olat), np.nanmax(olat))
                #: The real max LON/LAT can be obtain by
                #: np.max(lonremap_mesh[~np.isnan(ocaMap_ctp)])
                if no_clas_data:
                    claasMap_ctp = np.zeros(lonremap_mesh.shape) + np.nan
                    claasMap_cfc = np.zeros(lonremap_mesh.shape) + np.nan
                    claasExtent = (minLon, maxLon, minLat, maxLat)
                else:
                    claasMap_ctp = np.where(clasIndctp, cctp, np.nan)
                    claasMap_cfc = np.where(clasIndcfc, ccfc, np.nan)
                    claasExtent = (np.nanmin(ctolon), np.nanmax(ctolon), np.nanmin(ctolat), np.nanmax(ctolat))
                
                if no_gewex_data:
                    gewexMap_ctp = np.zeros(lonremap_mesh.shape) + np.nan
                    gewexMap_cfc = np.zeros(lonremap_mesh.shape) + np.nan
                    gewexExtent = (minLon, maxLon, minLat, maxLat)
                else:
                    gewexMap_ctp = np.where(gewexIndctp, gctp, np.nan)
                    gewexMap_cfc = np.where(gewexIndcfc, gcfc, np.nan)
                    gewexExtent = (np.nanmin(glon), np.nanmax(glon), np.nanmin(glat), np.nanmax(glat))
                plotMaps_done = True
#     class cartopy.crs.Orthographic(central_longitude=0.0, central_latitude=0.0, globe=None)[source]
# 
# 
# 
    
#     octp_use = []
#     cctp_use = []
#     for i in range(olon.shape[0]-1):
    if args.om:
        sys.exit()
    oca_ctp_plot = oca_ctp[0:i+1]
    clas_ctp_plot = clas_ctp[0:i+1]
    gewex_ctp_plot = gewex_ctp[0:i+1]
    oca_cfc_plot = oca_cfc[0:i+1]
    clas_cfc_plot = clas_cfc[0:i+1]
    gewex_cfc_plot = gewex_cfc[0:i+1]
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
    
    if use_opt in [12, 22]:
        title_extra = '\n Area = +-%d, Resample = True' %(int(np.nanmax(olon[ocaIndctp]) + 0.5))
        
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
        ax.set_xticks(plot_x_ticks)
        ax.set_xticklabels(plot_x_labels)
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_ylabel('Cloud Top Pressure [hPa]')
        ax.set_ylim([700, 400])
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
        ax.set_ylim([45, 75])
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
        im = ax.imshow(ocaMap_ctp, origin='lower', extent = ocaExtent, transform=ccrs.PlateCarree(), vmin=100, vmax=1000)#, alpha=0.5)
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
        im = ax.imshow(claasMap_ctp, origin='lower', extent = claasExtent, transform=ccrs.PlateCarree(), vmin=100, vmax=1000)#, alpha=0.5)
        ax.coastlines(resolution='110m')#, color='r')
        ax.set_title('CLAAS, %d-%02d' %(year, mon) + title_extra)
        cbar = fig.colorbar(im)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
        cbar.set_label('Cloud Top Pressure [hPa]')
        figname = '%s/map_ctp_claas_month-mean-%d-%02d' %(plotDir, year, mon) + '_opt-%d' %use_opt
        fig.savefig(figname +'.png')
        print(figname)
        
        #: Gewex is not available for all months
        if not np.isnan(gewexMap_ctp).all():
            plt.figure(figsize=(3, 3))
            fig = plt.figure()
            ax = fig.add_subplot(111,projection=ccrs.Orthographic())#central_longitude=-80))
        #     ax.gridlines()
            im = ax.imshow(gewexMap_ctp, origin='lower', extent = gewexExtent, transform=ccrs.PlateCarree(), vmin=100, vmax=1000)#, alpha=0.5)
        #     ax.imshow(g_ctp_ocaLL, origin='lower', extent = (maxLon, minLon, maxLat, minLat), transform=ccrs.PlateCarree())#, alpha=0.5)
            ax.coastlines(resolution='110m')#, color='r')
            ax.set_title('GEWEX, %d-%02d' %(year, mon) + title_extra)
            cbar = fig.colorbar(im)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
            cbar.set_label('Cloud Top Pressure [hPa]')
            figname = '%s/map_ctp_gewex_month-mean-%d-%02d' %(plotDir, year, mon) + '_opt-%d' %use_opt
            fig.savefig(figname +'.png')
            print(figname)
        
        
        
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
        
        
        #: CFC
        
        plt.figure(figsize=(3, 3))
        fig = plt.figure()
        ax = fig.add_subplot(111,projection=ccrs.Orthographic())
    #     ax.gridlines()
        im = ax.imshow(ocaMap_cfc, origin='lower', extent = ocaExtent, vmin=0, vmax=100, transform=ccrs.PlateCarree())#, alpha=0.5)
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
        im = ax.imshow(claasMap_cfc, origin='lower', extent = claasExtent, transform=ccrs.PlateCarree(), vmin=0, vmax=100)#, alpha=0.5)
        ax.coastlines(resolution='110m')#, color='r')
        ax.set_title('CLAAS, %d-%02d' %(year, mon) + title_extra)
        cbar = fig.colorbar(im)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
        cbar.set_label('Cloud Fraction [%]')
        figname = '%s/map_cfc_claas_month-mean-%d-%02d' %(plotDir, year, mon) + '_opt-%d' %use_opt
        fig.savefig(figname +'.png')
        print(figname)
        
        #: Gewex is not available for all months
        if not np.isnan(gewexMap_cfc).all():
            plt.figure(figsize=(3, 3))
            fig = plt.figure()
            ax = fig.add_subplot(111,projection=ccrs.Orthographic())#central_longitude=-80))
        #     ax.gridlines()
            im = ax.imshow(gewexMap_cfc, origin='lower', extent = gewexExtent, transform=ccrs.PlateCarree(), vmin=0, vmax=100)#, alpha=0.5)
        #     ax.imshow(g_ctp_ocaLL, origin='lower', extent = (maxLon, minLon, maxLat, minLat), transform=ccrs.PlateCarree())#, alpha=0.5)
            ax.coastlines(resolution='110m')#, color='r')
            ax.set_title('GEWEX, %d-%02d' %(year, mon) + title_extra)
            cbar = fig.colorbar(im)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
            cbar.set_label('Cloud Fraction [%]')
            figname = '%s/map_cfc_gewex_month-mean-%d-%02d' %(plotDir, year, mon) + '_opt-%d' %use_opt
            fig.savefig(figname +'.png')    
            print(figname)





    
    
    
    
    
    
    