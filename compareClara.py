#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Created on 2020-08-14

Copyright (c) 2020 Erik Johansson

Author(s):
    Erik Johansson <erik.johansson@smhi.se>
'''

import numpy as np
import netCDF4
import pdb
import glob
from matplotlib import pyplot as plt  # @UnresolvedImport
import datetime
import cartopy.crs as ccrs  # @UnresolvedImport
import sys

import matplotlib  # @UnresolvedImport
from pyresample import geometry as prg  # @UnresolvedImport
from pyresample.kd_tree import resample_nearest  # @UnresolvedImport
import copy
from scipy.ndimage.interpolation import zoom
import time

CB_colour_cycle = {'blue': '#377eb8', 'orange': '#ff7f00', 'green': '#4daf4a',
                   'pink': '#f781bf', 'brown': '#a65628', 'lila': '#984ea3',
                   'grey': '#999999', 'red': '#e41a1c', 'gul': '#dede00'}


#: Figure variables
def getFigVar(vn):
    if vn == 'cfc':
        mi_o = 0
        ma_o = 100
        bt_o = [0, 20, 40, 60, 80, 100]
        cl_o = '[%]'
        
        mi_d = -30
        ma_d = 30
        bt_d = [-30, -15, 0, 15, 30]
        cl_d = '[%%]'
        
    elif vn == 'cth':
        mi_o = 0
        ma_o = 16000
        bt_o = [2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000]
        ma_o = 10000
        bt_o = [2000, 4000, 6000, 8000, 10000]#, 12000, 14000, 16000]
        cl_o = '[m]'
        
        mi_d = -5000
        ma_d = 5000
        bt_d = [-5000, -2500, 0, 2500, 5000]
        cl_d = '[m]'
    elif vn == 'cph':
        mi_o = 0
        ma_o = 100
        bt_o = [0, 20, 40, 60, 80, 100]
        cl_o = '[%]'
        
        mi_d = -30
        ma_d = 30
        bt_d = [-30, -15, 0, 15, 30]
        cl_d = '[%%]'
    elif vn == 'lwp':
        mi_o = 0
        ma_o = 0.25
        bt_o = [0, 0.05, 0.1, 0.15, 0.2, 0.25]
        cl_o = '[kg/m2]'
        
        mi_d = -0.1
        ma_d = mi_d * -1
        bt_d = [-0.1, -0.05, 0, 0.05, 0.1]
        cl_d = '[kg/m2]'
    elif vn == 'iwp':
        mi_o = 0
        ma_o = 0.25
        bt_o = [0, 0.05, 0.1, 0.15, 0.2, 0.25]
        cl_o = '[kg/m2]'
        
        mi_d = -0.1
        ma_d = mi_d * -1
        bt_d = [-0.1, -0.05, 0, 0.05, 0.1]
        cl_d = '[kg/m2]'
    elif vn == 'jch':
        mi_o = 0
        ma_o = 100
        bt_o = [0, 20, 40, 60, 80, 100]
        cl_o = '[%]'
        
        mi_d = -20
        ma_d = 20
        bt_d = [-20, -10, 0, 10, 20]
        cl_d = '[%%]'
    else:
        mi_o = 0
        ma_o = 100
        bt_o = [0, 20, 40, 60, 80, 100]
        cl_o = '[%]'
        
        mi_d = -30
        ma_d = 30
        bt_d = [-30, -15, 0, 15, 30]
        cl_d = '[%%]'
        
    retv = {'vmino': mi_o, 'vmaxo': ma_o, 'bartickso': bt_o, 'cbarlabelo': cl_o, \
            'vmind': mi_d, 'vmaxd': ma_d, 'barticksd': bt_d, 'cbarlabeld': cl_d}
    return retv


def getVar(ncO, vn, geo=False, obtv=None):
    #: Use obtv if variable should be added to already created dictionary
    if obtv is None:
        retv = {}
    else:
        retv = obtv
    if geo:
        #: Geo location
        time = np.where(ncO.variables['time'][:].mask, np.nan, ncO.variables['time'][:].data)
        lat = np.where(ncO.variables['lat'][:].mask, np.nan, ncO.variables['lat'][:].data)
        lon = np.where(ncO.variables['lon'][:].mask, np.nan, ncO.variables['lon'][:].data)
        latm = ncO.variables['lat'][:]
        lonm = ncO.variables['lon'][:]
        retv.update({'time': time, 'lat': lat, 'lon': lon, 'latm': latm, 'lonm': lonm})
    if vn == 'jch':
        vnr = 'cfc'
    else:
        vnr = vn
    if vnr not in retv.keys():
        #: Do not add variables that are already added. Ex lat if geo is True
#         pdb.set_trace()
        v = np.where(ncO.variables[vnr][:].mask, np.nan, ncO.variables[vnr][:].data)
        vm = ncO.variables[vnr][:]
        retv.update({vn: v, (vn + 'm'): vm})
    return retv


if __name__ == '__main__':
#     mainStange = '/nobackup/smhid15/sm_kgkar/CLARA_A3_feedbackloop/strange_cases'
#     sfile = mainStange + '/' + 'S_NWC_CMAPROB_noaa15_54862_20081201T1048395Z_20081201T1238305Z.nc'
    mainClaraA3Dir = '/nobackup/smhid15/sm_kgkar/CLARA_A3_feedbackloop'
    mainClaraA2Dir = '/nobackup/smhid11/foua/data/CLARA-A2_final'
    mainClaraA21Dir = '/nobackup/smhid11/foua/data/CLARA-A2.1'
    mainCalipsoDir = '/nobackup/smhid17/proj/foua/data/satellit/calipso/Cloud_Occurrence_V1'
    
    plotDir = '/nobackup/smhid17/users/sm_erjoh/Cmsaf/Plots/Compare_A3_A2'
    
    yearmonths = [197901, 197907, 198104, 198301, 198307, 198701, 198707, 199201, 199207, \
                  199501, 199507, 199801, 199807, 200204, 200801, 200802, 200803, 200804, \
                  200805, 200806, 200807, 200808, 200809, 200810, 200811, 200812, 200904, \
                  201001, 201007, 201101, 201107, 201601, 201607]
    
    #: 'CFCmm', 'CTOmm', 'CPHmm', 'LWPmm', 'IWPmm', 'JCHmh'
    typ = 'CFCmm'
    var = typ[0:3].lower()
    if var == 'cto':
        var = 'cth'
    fv = getFigVar(var)
    for area in ['GL', 'NP', 'SP']:
        cfc3mean = []
        cfc2mean = []
        ym_d3 = []
        ym_d2 = []
        for ym in yearmonths:
            year = int(str(ym)[0:4])
#             if year != 2008:
#                 continue
            mon = int(str(ym)[4:])
            if (area in ['NP', 'SP']) and (typ in ['CTOmm', 'CPHmm', 'LWPmm', 'IWPmm', 'JCHmh']):
                continue
#             if (area in ['NP', 'SP']) and (year == 2016):
#                 continue
#             if (var == 'jch') and (year == 2016):
#                 continue
#             if year < 1982:
#                 continue
            if year in [2016]:
                c2d = '%s/*/AVPOS/%d' %(mainClaraA21Dir, year)
            else:
                c2d = '%s/%d/%02d/nc/AVPOS_%d_GAC_V002_L3' %(mainClaraA2Dir, year, mon, ym)
            print(ym)
            c3d = '%s/%d/AVPOS_%d_CLARA3_Level3_V011' %(mainClaraA3Dir, ym, ym)
            fn3 = glob.glob('%s/%s%d*%s.nc' %(c3d, typ, ym, area))
            fn2 = glob.glob('%s/%s%d*%s.nc' %(c2d, typ, ym, area))
            if (len(fn3) + len(fn2)) == 0:
                print('No files')
                continue
            if len(fn3) == 1:
                nc3 = netCDF4.Dataset(fn3[0])
                val3 = getVar(nc3, var, geo=True)
                nc3.close()
                ym_d3.append(datetime.datetime(year=year, month=mon, day=1))
                if (area == 'GL') and (var in ['iwp', 'lwp']):
                    lat3 = (np.abs(val3['lat']) <= 60)
                    cfc3mean.append(np.nanmean(val3[var][0,lat3,:]))
                else:
                    cfc3mean.append(np.nanmean(val3[var]))
            elif len(fn3) > 1:
                print('To many a3 files')
                pdb.set_trace()
            if len(fn2) == 1:
                nc2 = netCDF4.Dataset(fn2[0])
                val2 = getVar(nc2, var, geo=True)
                nc2.close()
                ym_d2.append(datetime.datetime(year=year, month=mon, day=1))
                if (area == 'GL') and (var in ['iwp', 'lwp']):
                    lat2 = (np.abs(val2['lat']) <= 60)
                    cfc2mean.append(np.nanmean(val2[var][0,lat2,:]))
                else:
                    cfc2mean.append(np.nanmean(val2[var]))
            elif len(fn2) > 1:
                print('To many a2 files')
                pdb.set_trace()
            
            
#             #: Calipso
#             if year == 2008:
#                 cald = '%s/%d/%02d' %(mainCalipsoDir, year, mon)
#                 calf = glob.glob('%s/CAL_LID_L3_Cloud_Occurrence-Standard-V1-00.%d-%02dA.h5' %(cald, year, mon))[0]
#                 import h5py
#                 h5cal = h5py.File(calf, 'r')
#             
#                 pdb.set_trace()
            if (ym in [198301, 201107]) or (year == 2008):
                #: The data is in lon/lat therefore use PlateCarree?
                img_proj = ccrs.PlateCarree()
                if area == 'GL':
                    #: Extent of the data
                    img_extent = (-180, 180, -90, 90)
                    zoom_val = 1
                    #: Projection of figure
                    projection = ccrs.PlateCarree()
                    data3 = val3[var][0,:,:]
                    data2 = val2[var][0,:,:]
                    #: lon/lat
                    lons, lats = np.meshgrid(val3['lon'], val3['lat'])
                    #: figsize
                    fs = (10,20)
                elif area in ['NP', 'SP']:
                    if area == 'NP':
                        img_extent = (-180, 180, 45, 90)
                        projection = ccrs.Orthographic(central_longitude=0.0, central_latitude=90.0)#, globe=None)
                    else:
                        img_extent = (-180, 180, -45, -90)
                        projection = ccrs.Orthographic(central_longitude=0.0, central_latitude=-90.0)#, globe=None)
                        
                    zoom_val = 0.75
#                     data3 = np.flipud(val3[var][0,:,:])
                    data3 = np.fliplr(val3[var][0,:,:])
#                     data3 = val3[var][0,:,:]
                    data2 = np.delete(val2[var][0,:,:], (180), axis=0)
                    data2 = np.delete(data2, (180), axis=1)
                    data2 = np.fliplr(data2)
                    lons = val3['lon']
                    lats = val3['lat']
                    fs = (5,20)

                datad = (data3 - data2)
                
                #: Create colormap with nan's not visible
                my_cmap = copy.copy(matplotlib.cm.cividis_r)
                my_cmap.set_bad('1.0', alpha=0)
#                 my_cmap.set_over('red', alpha=1)                

                #: Define area from cartopy to get globe boarders
                crs = projection
                crs.bounds = (crs.x_limits[0]*zoom_val, crs.x_limits[1]*zoom_val, \
                              crs.y_limits[0]*zoom_val, crs.y_limits[1]*zoom_val)
                
                #: Create pyresample area_def object for resampling
                area_def = prg.AreaDefinition('orthographic',
                                              'orthographic',
                                              'orthographic',
                                              projection=crs.proj4_params,
                                              width = lons.shape[1], height = lons.shape[0], 
                                              area_extent=(crs.x_limits[0]*zoom_val,
                                                           crs.y_limits[0]*zoom_val,
                                                           crs.x_limits[1]*zoom_val,
                                                           crs.y_limits[1]*zoom_val))
                
                
                # Remap to crs projection
                swath_def = prg.SwathDefinition(lons=lons, lats=lats)
                result3 = resample_nearest(swath_def, data3, area_def, 
                                           radius_of_influence=lons.shape[0] * lons.shape[1] * 2.5, fill_value=None)
                result2 = resample_nearest(swath_def, data2, area_def, 
                                           radius_of_influence=lons.shape[0] * lons.shape[1] * 2.5, fill_value=None)
                resultd = resample_nearest(swath_def, datad, area_def, 
                                           radius_of_influence=lons.shape[0] * lons.shape[1] * 2.5, fill_value=None)
                # Hack to remove data in corners
                lon, lat = area_def.get_lonlats()
                yshape, xshape = result3.data.shape
                result3.data[:,0:xshape//2][lon[:,0:xshape//2]>0] = np.nan
                result3.data[:,xshape//2:][lon[:,xshape//2:]<0] = np.nan
                
                result2.data[:,0:xshape//2][lon[:,0:xshape//2]>0] = np.nan
                result2.data[:,xshape//2:][lon[:,xshape//2:]<0] = np.nan
                
                resultd.data[:,0:xshape//2][lon[:,0:xshape//2]>0] = np.nan
                resultd.data[:,xshape//2:][lon[:,xshape//2:]<0] = np.nan
                
                fig = plt.figure(figsize=fs)
                fig.suptitle('%s %d-%02d' %(var.upper(), year, mon))
                ax = fig.add_subplot(3,1,1, projection=crs)
                im = ax.imshow(result3, transform=crs, extent=crs.bounds, cmap=my_cmap, vmin=fv['vmino'], vmax=fv['vmaxo'])
                ax.coastlines()
                cbar = fig.colorbar(im, ticks=fv['bartickso'])
                cbar.set_label(fv['cbarlabelo'], rotation=90)
                ax.set_title('CLARA3')
                
                ax = fig.add_subplot(3,1,2, projection=crs)
                im = ax.imshow(result2, transform=crs, cmap=my_cmap, extent=crs.bounds, vmin=fv['vmino'], vmax=fv['vmaxo'])
                ax.coastlines()
                cbar = fig.colorbar(im, ticks=fv['bartickso'])
                cbar.set_label(fv['cbarlabelo'], rotation=90)
                ax.set_title('CLARA2')
                
#                 vminmaxd = 90
#                 barticksd = [-90, -60, -30, 0, 30, 60, 90]
                ax = fig.add_subplot(3,1,3, projection=crs)
                im = ax.imshow(resultd, transform=crs, extent=crs.bounds, cmap='RdBu_r', vmin=fv['vmind'], vmax=fv['vmaxd'])
                ax.coastlines()
                cbar = fig.colorbar(im, ticks=fv['barticksd'])
                cbar.set_label(fv['cbarlabeld'], rotation=90)
                ax.set_title('A3 - A2')
                
                figname = '%s/%s_%s_map_%d-%02d' %(plotDir, var, area, year,mon)
                fig.show()
                pdb.set_trace()
                fig.savefig(figname + '.png')

        
        if len(cfc3mean) == 0:
            continue
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(ym_d3, cfc3mean, color=CB_colour_cycle['red'], label='A3')
        ax.plot(ym_d2, cfc2mean, color=CB_colour_cycle['blue'], label='A2')
        #: Check where both has value
        d_3 = []
        d_2 = []
        for i3, val in enumerate(ym_d3):
            if val in ym_d2:
                i2 = np.where(np.asarray(ym_d2) == val)[0][0]
                d_3.append(cfc3mean[i3])
                d_2.append(cfc2mean[i2])
        #: Calculate difference
        d_u = np.asarray(d_3) - np.asarray(d_2)
        d_um = np.mean(d_u)
        d_p = (d_u / np.asarray(d_2)) * 100
        d_pm = np.mean(d_p)
    #     for i, txt in enumerate(proc):
#     #         ax.annotate(txt, (cfc3Gmean[i], ym_d[i]))
     
     
     
        fig.autofmt_xdate()
        ax.set_xlabel('Date')
        if var in ['iwp', 'lwp']:
            ax.set_ylabel('%s Mean (lat $\pm$ 60) %s %s' %(area, var.upper(), fv['cbarlabelo']))
            ax.text(0.5, 0.9, 'Mean increase = %0.3f%s, %0.1f%%' %(d_um, fv['cbarlabeld'].replace('[', '').replace(']', ''), d_pm), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        else:
            ax.set_ylabel('%s Mean %s %s' %(area, var.upper(), fv['cbarlabelo']))
            ax.text(0.5, 0.9, 'Mean increase = %0.1f%s, %0.1f%%' %(d_um, fv['cbarlabeld'].replace('[', '').replace(']', ''), d_pm), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        ax.legend()
        fig.tight_layout()
        figname = '%s/%s_%s_mean' %(plotDir, var, area)
        fig.show()
        fig.savefig(figname + '.png')
    pdb.set_trace()
    
    
    
    
    
    
    
    
    
    
    
    