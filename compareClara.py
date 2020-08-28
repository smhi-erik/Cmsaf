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

CB_colour_cycle = {'blue': '#377eb8', 'orange': '#ff7f00', 'green': '#4daf4a',
                   'pink': '#f781bf', 'brown': '#a65628', 'lila': '#984ea3',
                   'grey': '#999999', 'red': '#e41a1c', 'gul': '#dede00'}


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
    if vn not in retv.keys():
        #: Do not add variables that are already in added. Ex lat if geo is True
        v = np.where(ncO.variables[vn][:].mask, np.nan, ncO.variables[vn][:].data)
        vm = ncO.variables[vn][:]
        retv.update({vn: v, (vn + 'm'): vm})
    return retv


if __name__ == '__main__':
#     mainStange = '/nobackup/smhid15/sm_kgkar/CLARA_A3_feedbackloop/strange_cases'
#     sfile = mainStange + '/' + 'S_NWC_CMAPROB_noaa15_54862_20081201T1048395Z_20081201T1238305Z.nc'
    mainClaraA3Dir = '/nobackup/smhid15/sm_kgkar/CLARA_A3_feedbackloop'
    mainClaraA2Dir = '/nobackup/smhid11/foua/data/CLARA-A2_final'
    mainClaraA21Dir = '/nobackup/smhid11/foua/data/CLARA-A2.1'
    
    plotDir = '/nobackup/smhid17/users/sm_erjoh/Cmsaf/Plots/Compare_A3_A2'
    
    yearmonths = [197901, 197907, 198104, 198301, 198307, 198701, 198707, 199201, 199207, \
                  199501, 199507, 199801, 199807, 200204, 200801, 200802, 200803, 200804, \
                  200805, 200806, 200807, 200808, 200809, 200810, 200811, 200812, 200904, \
                  201001, 201007, 201101, 201107, 201601, 201607]
    
    #: 'CFCmm', 'CTOmm', 'CPHmm', 'LWPmm', 'IWPmm', 'JCHmh'
    typ = 'CFCmm'
    var = typ[0:3].lower()
    for area in ['GL', 'NP', 'SP']:
        cfc3mean = []
        cfc2mean = []
        ym_d = []
        for ym in yearmonths:
            year = int(str(ym)[0:4])
            mon = int(str(ym)[4:])
            if (area in ['NP', 'SP']) and (typ in ['CPHmm', 'LWPmm', 'IWPmm', 'JCHmh']):
                continue
            if (area in ['NP', 'SP']) and (year == 2016):
                continue
            if year < 1982:
                continue
            elif year in [2016]:
                c2d = '%s/*/AVPOS/%d' %(mainClaraA21Dir, year)
            else:
                c2d = '%s/%d/%02d/nc/AVPOS_%d_GAC_V002_L3' %(mainClaraA2Dir, year, mon, ym)
            print(ym)
            ym_d.append(datetime.datetime(year=year, month=mon, day=1))
        
            c3d = '%s/%d/AVPOS_%d_CLARA3_Level3_V011' %(mainClaraA3Dir, ym, ym)
        
            fn3 = glob.glob('%s/%s%d*%s.nc' %(c3d, typ, ym, area))[0]
            fn2 = glob.glob('%s/%s%d*%s.nc' %(c2d, typ, ym, area))[0]
            
            nc3 = netCDF4.Dataset(fn3)
            nc2 = netCDF4.Dataset(fn2)
            val3 = getVar(nc3, var, geo=True)
            val2 = getVar(nc2, var, geo=True)
            nc3.close()
            nc2.close()

            
            
            cfc3mean.append(np.nanmean(val3[var]))
            cfc2mean.append(np.nanmean(val2[var]))
            
            if (ym in [198301, 201107]):
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
                my_cmap = copy.copy(matplotlib.cm.BrBG)
                my_cmap.set_bad('1.0', alpha=0)
                my_cmap.set_over('red', alpha=1)
                

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
                im = ax.imshow(result3, transform=crs, extent=crs.bounds, cmap=my_cmap, vmin=0, vmax=100)
                ax.coastlines()
                barticks = [0, 20, 40, 60, 80, 100]
                cbar = fig.colorbar(im, ticks=barticks)
                cbar.set_label('[%]', rotation=0)
                ax.set_title('CLARA3')
                
                ax = fig.add_subplot(3,1,2, projection=crs)
                im = ax.imshow(result2, transform=crs, cmap=my_cmap, extent=crs.bounds, vmin=0, vmax=100)
                ax.coastlines()
                barticks = [0, 20, 40, 60, 80, 100]
                cbar = fig.colorbar(im, ticks=barticks)
                cbar.set_label('[%]', rotation=0)
                ax.set_title('CLARA2')
                
                ax = fig.add_subplot(3,1,3, projection=crs)
                im = ax.imshow(resultd, transform=crs, extent=crs.bounds, cmap='RdBu_r', vmin=-90, vmax=90)
                ax.coastlines()
                barticks = [-90, -60, -30, 0, 30, 60, 90]
                cbar = fig.colorbar(im, ticks=barticks)
                ax.set_title('A3 - A2')
                
                figname = '%s/%s_%s_map_%d-%02d' %(plotDir, var, area, year,mon)
                fig.show()
                pdb.set_trace()
                fig.savefig(figname + '.png')
                

                
#                 fig = plt.figure()
#                 ax = fig.add_subplot(1,1,1, projection=projection)
# #                 fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection=projection))
#                 ax.set_extent([-180, 180, 45, 90], ccrs.PlateCarree())
#                 im = ax.pcolormesh(val3['lonm'].data, val3['latm'].data, np.fliplr(val3['cfcm'].data[0,:,:]), vmin=0, vmax=100.0, cmap=cmap, transform=ccrs.PlateCarree())
# #                 contour_levels = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.]
# #                 im = ax.contourf(val3['lonm'].data, val3['latm'].data, val3['cfcm'].data[0,:,:], \
# #                                  vmin=0, vmax=100.0, cmap=cmap, transform=ccrs.PlateCarree(), levels=contour_levels)
# #                                  levels=contour_levels, transform=img_proj, extent=img_extent)
# #                 barticks = [0, 20, 40, 60, 80, 100]
# #                 cbar = fig.colorbar(im, ticks=barticks)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
#                 cbar = fig.colorbar(im)
#                 ax.coastlines()
#                 fig.savefig('%s/%s_%s_map_pcm.png' %(plotDir, var, area))
#                 fig.show()
#                 pdb.set_trace()
#                 sys.exit()
                
#                 figa = plt.figure(figsize=(10,20))
#                 figa.suptitle('pcolormesh')#'%s %d-%02d' %(var.upper(), year, mon))
#                 axa1 = figa.add_subplot(3,1,1, projection=projection)
#                 axa1.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
#                 axa1.coastlines()
# #                 axa1.gridlines()
#                 ima1 = axa1.pcolormesh(val3['lon'], val3['lat'], data3, cmap=my_cmap, vmin=0, vmax=100, \
#                                      transform=img_proj)
#                 barticks = [0, 20, 40, 60, 80, 100]
#                 cbara1 = figa.colorbar(ima1, ticks=barticks)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
#                 cbara1.set_label('[%]', rotation=0)
#                 axa1.set_title('CLARA3')
#                 
#                 axa2 = figa.add_subplot(3,1,2, projection=ccrs.PlateCarree())
#                 axa2.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
#                 axa2.coastlines()
#                 ima2 = axa2.pcolormesh(val2['lon'], val2['lat'], data2, cmap=my_cmap, vmin=0, vmax=100, \
#                                      transform=img_proj)
#                 cbara2 = figa.colorbar(ima2, ticks=barticks)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
#                 cbara2.set_label('[%]', rotation=0)
#                 axa2.set_title('CLARA2')
#                 
#                 axa3 = figa.add_subplot(3,1,3, projection=ccrs.PlateCarree())
#                 axa3.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
#                 axa3.coastlines()
#                 ima3 = axa3.pcolormesh(val3['lon'], val3['lat'], datad, cmap='RdBu_r', vmin=-90, vmax=90, \
#                                      transform=img_proj)
#                 barticks = [-90, -60, -30, 0, 30, 60, 90]
#                 cbara3 = figa.colorbar(ima3, ticks=barticks)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
#                 axa3.set_title('A3 - A2')
# #                 figa.tight_layout()
# #                 figname = '%s/%s_%s_map_%d-%02d' %(plotDir, var, area, year,mon)
#                 figa.show()
#                 pdb.set_trace()
#                 figa.savefig('%s/pcolormesh.png' %(plotDir))
#                 
#                 
#                 figa = plt.figure(figsize=(10,20))
#                 figa.suptitle('contourf')#'%s %d-%02d' %(var.upper(), year, mon))
#                 axa1 = figa.add_subplot(3,1,1, projection=projection)
#                 axa1.coastlines()
# #                 axa1.gridlines()
#                 contour_levels = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.]
#                 ima1 = axa1.contourf(val3['lon'], val3['lat'], val3[var][0,:,:], cmap=my_cmap, vmin=0, vmax=100, \
#                                      levels=contour_levels, transform=img_proj, extent=img_extent)
#                 barticks = [0, 20, 40, 60, 80, 100]
#                 cbara1 = figa.colorbar(ima1, ticks=barticks)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
#                 cbara1.set_label('[%]', rotation=0)
#                 axa1.set_title('CLARA3')
#                 
#                 axa2 = figa.add_subplot(3,1,2, projection=ccrs.PlateCarree())
#                 axa2.coastlines()
#                 ima2 = axa2.contourf(val2['lon'], val2['lat'], val2[var][0,:,:], cmap=my_cmap, vmin=0, vmax=100, \
#                                      levels=contour_levels, transform=img_proj, extent=img_extent)
#                 cbara2 = figa.colorbar(ima2, ticks=barticks)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
#                 cbara2.set_label('[%]', rotation=0)
#                 axa2.set_title('CLARA2')
#                 
#                 axa3 = figa.add_subplot(3,1,3, projection=ccrs.PlateCarree())
#                 axa3.coastlines()
#                 contour_levels = [-90., -75., -60., -45., -30., -15.,   0.,  15.,  30.,  45.,  60., 75., 90.]
#                 ima3 = axa3.contourf(val3['lon'], val3['lat'], datad, cmap='RdBu_r', vmin=-90, vmax=90, \
#                                      levels=contour_levels, transform=img_proj, extent=img_extent)
#                 cbara3 = figa.colorbar(ima3)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)
#                 axa3.set_title('A3 - A2')
# #                 figa.tight_layout()
# #                 figname = '%s/%s_%s_map_%d-%02d' %(plotDir, var, area, year,mon)
#                 figa.show()
#                 pdb.set_trace()
#                 figa.savefig('%s/contourf.png' %(plotDir))
# #                 figa.savefig(figname + '.png')
#                 pdb.set_trace()
    
    
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(ym_d, cfc3mean, color=CB_colour_cycle['red'], label='A3')
        ax.plot(ym_d, cfc2mean, color=CB_colour_cycle['blue'], label='A2')
        proc = ((np.asarray(cfc3mean) - np.asarray(cfc2mean))/np.asarray(cfc2mean)) * 100
        procmean = np.mean(proc)
    #     for i, txt in enumerate(proc):
    #         ax.annotate(txt, (cfc3Gmean[i], ym_d[i]))
        ax.text(0.5, 0.9, 'Mean increase = %0.2f%%' %procmean, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    
    
    
        fig.autofmt_xdate()
        ax.set_xlabel('Date')
        ax.set_ylabel('%s Mean %s' %(area, var.upper()))
        ax.legend()
        fig.tight_layout()
        figname = '%s/%s_%s_mean' %(plotDir, var, area)
        fig.show()
        fig.savefig(figname + '.png')  
        pdb.set_trace()
    
    
    
    
    
    
    
    
    
    
    
    