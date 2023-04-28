#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Created on 2023-03-30

Copyright (c) 2023 Erik Johansson

@author:     Erik Johansson
@contact:    <erik.johansson@smhi.se>
 
'''

import numpy as np
import netCDF4
import pdb


from pyresample.geometry import SwathDefinition
from pyresample.kd_tree import resample_nearest
from pyresample.kd_tree import get_neighbour_info
from pyresample.kd_tree import get_sample_from_neighbour_info

#: https://stackoverflow.com/questions/67779374/deprecationwarning-np-bool
from warnings import filterwarnings
import glob
import datetime
import matplotlib.pyplot as plt
import time
from scipy.stats import gaussian_kde
import scipy
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')
# filterwarnings(action='ignore', category=DeprecationWarning, message='tostring() is deprecated.')



def get_sunz_correction(scene):#, REFL_BANDS):
    #: Modifyed apply_sunz_correction 
    #: Copyed from level1c4pps 
    #: /Erik
    """Apply sun zenith angle correciton to visual channels."""
    sza = scene['sunzenith']
    mu0 = np.cos(np.radians(sza))
    scaler = 24.35 / (2 * mu0 + np.sqrt(498.5225 * mu0 * mu0 + 1))
#     for band in REFL_BANDS:
#         if band not in scene:
#             continue
#         if scene[band].attrs['sun_zenith_angle_correction_applied'] == 'False':
#             scene[band].values = scene[band].values * scaler
#             scene[band].attrs['sun_zenith_angle_correction_applied'] = 'True'
    if scaler.ndim == 3:
        scaler = scaler[0,:,:]
    return scaler

def getChannel(sati, chn, ang_masked):
    no_result = True
    for n in [0, 1, 2, 3, 4, 5, 6, 7, 20, 21, 22, 23, 24, 25, 26]:
        imageN = 'image%d' %n
        if imageN in sati.variables.keys():
            if chn == sati[imageN].id_tag:
                no_result = False
                ret = sati[imageN][0,:,:]
                #: only for visable channels
                if 'ch_r' in sati[imageN].id_tag:
                    if sati[imageN].sun_zenith_angle_correction_applied == 'False':
#                         from level1c4pps import apply_sunz_correction
                        scaler = get_sunz_correction(sati)
                        ret = ret * scaler
                if ret.mask.ndim != 0:
                    ret.mask[ang_masked] = True
                else:
                    ret.mask = ang_masked
                break
    if no_result:
        print('No result for %s' %chn)
        return None
    else:
        return ret


def getTimePerScanline(ncf, st):
    secFromStart = np.linspace(ncf['time_bnds'][:][0][0]*24*3600, ncf['time_bnds'][:][0][1]*24*3600, num=ncf['lon'].shape[0])
    start_time = datetime.datetime(int(st[0:4]), int(st[4:6]), int(st[6:8]), int(st[9:11]), int(st[11:13]), int(st[13:15]), int(st[15]))
    scanLineTime = [start_time + datetime.timedelta(seconds=x) for x in secFromStart]
    return np.asarray(scanLineTime)


def findEdges(n19lat, n19lon, n19time, npplat, npplon, npptime, timeDiffMin):
    maxDiff = timeDiffMin * 60
    minTime = np.max([n19time[0], npptime[0]]) - datetime.timedelta(seconds=maxDiff)
    maxTime = np.min([n19time[-1], npptime[-1]]) + datetime.timedelta(seconds=maxDiff)
    
    n19Cut = (n19time >= minTime) & (n19time <= maxTime)
    nppCut = (npptime >= minTime) & (npptime <= maxTime)
    
#             n19Cut = np.zeros(n19lat.shape).astype(bool)
#             nppCut = np.zeros(npplat.shape).astype(bool)
#             
#             nppStart = np.argmin(np.sqrt((npplat - n19lat[0])**2))
#             nppCut[nppStart:] = True
#             
#             n19Cut[0:nppCut.sum()] = True
    
    return n19Cut, nppCut



def cutEdges(obts, edg):
    retv = []
    for obt in obts:
        if obt is None:
            retv.append(obt)
        elif obt.ndim == 1:
            retv.append(obt[edg])
        elif obt.ndim == 2:
            retv.append(obt[edg, :])
        else:
            print("Strange dim in cut")
            pdb.set_trace()
    return retv


def plotScatterHisto(x_flat, y_flat, stitle, satn, xylim, pt_int, pt_str, pt_loc, histd2_bins, vmax, figname):
    fig = plt.figure()
    fig.suptitle(stitle)
    ax = fig.add_subplot(1,2,1)
    ax.scatter(x_flat, y_flat)
    ax.set_title('Scatter')
    ax.set_xlabel('NOAA 19, max = %d' %int(x_flat.max()))
    ax.set_ylabel('%s, max = %d' %(satn, int(y_flat.max())))
    ax.set_xlim(xylim)
    ax.set_ylim(xylim)
    ax.set_xticks(pt_int)
    ax.set_yticks(pt_int)
    
    xn = np.linspace(xylim[0], xylim[1], 100)
    k2, m = np.ma.polyfit(x_flat, y_flat, 1)
    ny_y2 = k2*xn + m
    ax.plot(xn, ny_y2, 'g--', label = '%.4G*x%+.2f' %(k2, m))
    rr = np.ma.corrcoef(x_flat, y_flat)
    
#         ax.plot(xn, xn, 'r', label = '1*x')
    ax.text(1,0,'ccof = %.4f' %rr[0,1], horizontalalignment='right', transform=ax.transAxes)
    ax.legend(loc=2)
    ax.set_aspect(1)
    
    
    ax = fig.add_subplot(1,2,2)
    H, xedges, yedges = np.histogram2d(x_flat, y_flat, bins=histd2_bins, range=[xylim, xylim])
    cmap = 'plasma'
    im = ax.imshow(H.T, origin='lower', vmin=0, vmax=vmax, cmap=cmap)
    #: -1 since the lat one is included in contrast to range. -0.5 to get everything inside the plot
    xnH = np.linspace(0, histd2_bins-1, 100)
    ny_y2H = k2*xnH + m
    ny_y2H = np.where(ny_y2H >= 99.75, np.nan, ny_y2H)
    
    zi = (ny_y2 - xylim[0]) / (xylim[1] - xylim[0]) * 100
    zi = np.where(zi >= 99.75, np.nan, zi)
    ax.plot(xnH, zi, 'g--')
#         ax.plot(xnH, ny_y2H, 'g--')
    
    
#         ax.plot(xnH, xnH, 'r')
    ax.set_title('2D Histogram')
    ax.set_xticks(pt_loc)
    ax.set_xticklabels(pt_str)
    ax.set_yticks(pt_loc)
    ax.set_yticklabels(pt_str)
    fig.subplots_adjust(right=0.89)
    pos2 = ax.get_position()
    cbar_ax = fig.add_axes([0.90, pos2.y0, 0.01, pos2.y1 - pos2.y0])
    cbar = fig.colorbar(im, cax=cbar_ax)
#         fig.savefig(figname + '.png')
    fig.show()
    

if __name__ == '__main__':
    mainDir = 'Data/testdata1_VGAC_CLARA35_corrected'
    radianceDir = '%s/radiances_l1c' %mainDir
    plotDir = 'Plots/Compare_VGAC'
#     n19f = '%s/S_NWC_avhrr_noaa19_99999_20130103T0656486Z_20130103T0838481Z.nc' %radianceDir
#     vgacf = '%s/S_NWC_avhrr_vgacsnpp_00000_20130103T0635000Z_20130103T0817290Z.nc' %radianceDir
#     viirsf = '%s/S_NWC_viirs_npp_00000_20130103T0635000Z_20130103T0817290Z.nc' %radianceDir
    
    n19_StartTimes = ['20130103T0656486', '20130116T0611437', '20130116T0759087', '20130116T0935387']
#     npp_StartTimes = ['20130103T0635000', '20130116T0550000', '20130116T0732000', '20130116T0913000']
    npp_StartTimes = ['20130103T0635580', '20130116T0551290', '20130116T0732590', '20130116T0914290']
#     n19 = Dataset(n19f, mode='r')
    
    #: 180 = inf for the angles
    accept_satz_max = 30
    accept_sunz_max = 80
    accept_time_diff = 1 #: in minutes
    
    a = -1
    for i in [1]:#range(len(n19_StartTimes)):
        a = a + 1
        n19f = glob.glob('%s/S_NWC_avhrr_noaa19_*%s*.nc' %(radianceDir, n19_StartTimes[i]))[0]
        vgacf = glob.glob('%s/S_NWC_avhrr_vgacsnpp_*%s*.nc' %(radianceDir, npp_StartTimes[i]))[0]
        viirsf = glob.glob('%s/S_NWC_viirs_npp_*%s*.nc' %(radianceDir, npp_StartTimes[i]))[0]
        
        n19 = netCDF4.Dataset(n19f, 'r', format='NETCDF4')
        vgac = netCDF4.Dataset(vgacf, 'r', format='NETCDF4')
        viirs = netCDF4.Dataset(viirsf, 'r', format='NETCDF4')
        
        
#     ct = pps_nc1['cma'][:]
#         if accept_satz_max == 180:
        if accept_satz_max == 180:
            n19_satza_masked = np.zeros(n19['satzenith'][0,:,:].data.shape).astype(bool)
            vgac_satza_masked = np.zeros(vgac['satzenith'][0,:,:].data.shape).astype(bool)
            viirs_satza_masked = np.zeros(viirs['satzenith'][0,:,:].data.shape).astype(bool)
        else:
            n19_satza_masked = n19['satzenith'][0,:,:].data > accept_satz_max
            vgac_satza_masked = vgac['satzenith'][0,:,:].data > accept_satz_max
            viirs_satza_masked = viirs['satzenith'][0,:,:].data > accept_satz_max
        if accept_sunz_max == 180:
            n19_sunza_masked = np.zeros(n19['sunzenith'][0,:,:].data.shape).astype(bool)
            vgac_sunza_masked = np.zeros(vgac['sunzenith'][0,:,:].data.shape).astype(bool)
            viirs_sunza_masked = np.zeros(viirs['sunzenith'][0,:,:].data.shape).astype(bool)
        else:
            n19_sunza_masked = n19['sunzenith'][0,:,:].data > accept_sunz_max
            vgac_sunza_masked = vgac['sunzenith'][0,:,:].data > accept_sunz_max
            viirs_sunza_masked = viirs['sunzenith'][0,:,:].data > accept_sunz_max
        n19_masked = n19_satza_masked | n19_sunza_masked
        vgac_masked = vgac_satza_masked | vgac_sunza_masked
        viirs_masked = viirs_satza_masked | viirs_sunza_masked

        n19_lats =  n19['lat'][:]
        n19_lons =  n19['lon'][:]
        n19_r06 = getChannel(n19, 'ch_r06', n19_masked)
        n19_r09 = getChannel(n19, 'ch_r09', n19_masked)
        n19_r16 = getChannel(n19, 'ch_r16', n19_masked)
        n19_tb11 = getChannel(n19, 'ch_tb11', n19_masked)
        n19_tb12 = getChannel(n19, 'ch_tb12', n19_masked)
        n19_tb37 = getChannel(n19, 'ch_tb37', n19_masked)
        
        vgac_lats =  vgac['lat'][:]
        vgac_lons =  vgac['lon'][:]
        vgac_r06 = getChannel(vgac, 'ch_r06', vgac_masked)
        vgac_r09 = getChannel(vgac, 'ch_r09', vgac_masked)
        vgac_r16 = getChannel(vgac, 'ch_r16', vgac_masked)
        vgac_tb11 = getChannel(vgac, 'ch_tb11', vgac_masked)
        vgac_tb12 = getChannel(vgac, 'ch_tb12', vgac_masked)
        vgac_tb37 = getChannel(vgac, 'ch_tb37', vgac_masked)
        
        viirs_lats =  viirs['lat'][:]
        viirs_lons =  viirs['lon'][:]
        viirs_r06 = getChannel(viirs, 'ch_r06', viirs_masked)
        viirs_r09 = getChannel(viirs, 'ch_r09', viirs_masked)
        viirs_r16 = getChannel(viirs, 'ch_r16', viirs_masked)
        viirs_tb11 = getChannel(viirs, 'ch_tb11', viirs_masked)
        viirs_tb12 = getChannel(viirs, 'ch_tb12', viirs_masked)
        viirs_tb37 = getChannel(viirs, 'ch_tb37', viirs_masked)
        
        n19_scanLineTime = getTimePerScanline(n19, n19_StartTimes[i])
        npp_scanLineTime = getTimePerScanline(vgac, npp_StartTimes[i])
        
        npp_center_scanline = int(vgac_lats.shape[1]/2) + 1
        n19_center_scanline = int(n19_lats.shape[1]/2) + 1

        n19.close()
        vgac.close()
        viirs.close()
        
        if False:#i == 1:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(npp_scanLineTime, vgac_lats[:, npp_center_scanline])
            ax.plot(n19_scanLineTime, n19_lats[:, n19_center_scanline])
            fig.show()
        n19_use, npp_use = findEdges(n19_lats[:, n19_center_scanline], n19_lons[:, n19_center_scanline], n19_scanLineTime, vgac_lats[:, npp_center_scanline], vgac_lons[:, npp_center_scanline], npp_scanLineTime, accept_time_diff)
        
        n19_lats, n19_lons, n19_r06, n19_r09, n19_r16, n19_tb11, n19_tb12, n19_tb37 = cutEdges([n19_lats, n19_lons, n19_r06, n19_r09, n19_r16, n19_tb11, n19_tb12, n19_tb37], n19_use)
        vgac_lats, vgac_lons, vgac_r06, vgac_r09, vgac_r16, vgac_tb11, vgac_tb12, vgac_tb37 = cutEdges([vgac_lats, vgac_lons, vgac_r06, vgac_r09, vgac_r16, vgac_tb11, vgac_tb12, vgac_tb37], npp_use)
        viirs_lats, viirs_lons, viirs_r06, viirs_r09, viirs_r16, viirs_tb11, viirs_tb12, viirs_tb37 = cutEdges([viirs_lats, viirs_lons, viirs_r06, viirs_r09, viirs_r16, viirs_tb11, viirs_tb12, viirs_tb37], npp_use)
        
        if False:#i == 1:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(npp_scanLineTime[npp_use], vgac_lats[:, int(vgac_lats.shape[1]/2) + 1])
            ax.plot(n19_scanLineTime[n19_use], n19_lats[:, int(n19_lats.shape[1]/2) + 1])
            fig.show()
        
        
    
#     n19.close()
#     vgac.close()
        n19_swath_def = SwathDefinition(n19_lons, n19_lats)
        vgac_swath_def = SwathDefinition(vgac_lons, vgac_lats)
    #     euron1 = load_area("/home/a001865/git/acpg_develop/acpg/cfg/region_config.yaml", 'euron1')
        rm_n19_r06 = resample_nearest(n19_swath_def, n19_r06, vgac_swath_def,
                                  radius_of_influence=20000, fill_value=None)
        rm_n19_r09 = resample_nearest(n19_swath_def, n19_r09, vgac_swath_def,
                                  radius_of_influence=20000, fill_value=None)
        rm_n19_tb11 = resample_nearest(n19_swath_def, n19_tb11, vgac_swath_def,
                                  radius_of_influence=20000, fill_value=None)
        rm_n19_tb12 = resample_nearest(n19_swath_def, n19_tb12, vgac_swath_def,
                                  radius_of_influence=20000, fill_value=None)
        rm_n19_tb37 = resample_nearest(n19_swath_def, n19_tb37, vgac_swath_def,
                                  radius_of_influence=20000, fill_value=None)
#         
        
        
        if a == 0:
            rm_n19_r06_all = rm_n19_r06.copy()
            vgac_r06_all = vgac_r06.copy()
            viirs_r06_all = viirs_r06.copy()
            
            rm_n19_r09_all = rm_n19_r09.copy()
            vgac_r09_all = vgac_r09.copy()
            viirs_r09_all = viirs_r09.copy()
            
            rm_n19_tb11_all = rm_n19_tb11.copy()
            vgac_tb11_all = vgac_tb11.copy()
            viirs_tb11_all = viirs_tb11.copy()
            
            rm_n19_tb12_all = rm_n19_tb12.copy()
            vgac_tb12_all = vgac_tb12.copy()
            viirs_tb12_all = viirs_tb12.copy()
            
            rm_n19_tb37_all = rm_n19_tb37.copy()
            vgac_tb37_all = vgac_tb37.copy()
            viirs_tb37_all = viirs_tb37.copy()
        else:
            rm_n19_r06_all = np.ma.concatenate((rm_n19_r06_all, rm_n19_r06), axis=0)
            vgac_r06_all = np.ma.concatenate((vgac_r06_all, vgac_r06), axis=0)
            viirs_r06_all = np.ma.concatenate((viirs_r06_all, viirs_r06), axis=0)
            
            rm_n19_r09_all = np.ma.concatenate((rm_n19_r09_all, rm_n19_r09), axis=0)
            vgac_r09_all = np.ma.concatenate((vgac_r09_all, vgac_r09), axis=0)
            viirs_r09_all = np.ma.concatenate((viirs_r09_all, viirs_r09), axis=0)
            
            rm_n19_tb11_all = np.ma.concatenate((rm_n19_tb11_all, rm_n19_tb11), axis=0)
            vgac_tb11_all = np.ma.concatenate((vgac_tb11_all, vgac_tb11), axis=0)
            viirs_tb11_all = np.ma.concatenate((viirs_tb11_all, viirs_tb11), axis=0)
            
            rm_n19_tb12_all = np.ma.concatenate((rm_n19_tb12_all, rm_n19_tb12), axis=0)
            vgac_tb12_all = np.ma.concatenate((vgac_tb12_all, vgac_tb12), axis=0)
            viirs_tb12_all = np.ma.concatenate((viirs_tb12_all, viirs_tb12), axis=0)
            
            rm_n19_tb37_all = np.ma.concatenate((rm_n19_tb37_all, rm_n19_tb37), axis=0)
            vgac_tb37_all = np.ma.concatenate((vgac_tb37_all, vgac_tb37), axis=0)
            viirs_tb37_all = np.ma.concatenate((viirs_tb37_all, viirs_tb37), axis=0)
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     ax.scatter(rm_n19_r06_all, vgac_r06_all)
# #         ax.plot(n19_scanLineTime[n19_use], n19_lats[:, int(n19_lats.shape[1]/2) + 1])
#     fig.show()
    
    plotType = 'scatter'
#     plotType = 'density'
    r_max_axis = 180
    tb_min_axis = 170
    tb_max_axis = 350
    vmax = 100
    
    r_plot_ticks_int = [*range(0, r_max_axis+1, 30)]
    r_plot_ticks_str = np.asarray(r_plot_ticks_int).astype(str).tolist()
    r_plot_ticks_loc = []
    tb_plot_ticks_int = [*range(tb_min_axis, tb_max_axis+1, 30)]
    tb_plot_ticks_str = np.asarray(tb_plot_ticks_int).astype(str).tolist()
    tb_plot_ticks_loc = []
    histd2_bins = 100
    for i in range(len(r_plot_ticks_int)):
        #: -1 because we want the in between values
        r_plot_ticks_loc.append(int(np.round(i/(len(r_plot_ticks_int)-1) * histd2_bins)))
    for i in range(len(tb_plot_ticks_int)):
        tb_plot_ticks_loc.append(int(np.round(i/(len(tb_plot_ticks_int)-1) * histd2_bins)))
    #: Change last place to histd2_bins -1 instead of histd2_bins
    #: Othervise the figure looks vierd
    r_plot_ticks_loc[-1] = histd2_bins - 1
    tb_plot_ticks_loc[-1] = histd2_bins - 1
    
    title_end = ', SATZ < %d, SUNZ < %d, TD = %d min' %(accept_satz_max, accept_sunz_max, accept_time_diff)
    if accept_satz_max == 180:
        title_end = title_end.replace(', SATZ < 180', ', SATZ < inf')
    if accept_sunz_max == 180:
        title_end = title_end.replace(', SUNZ < 180', ', SATZ < inf')

    #: Channel r06
    rm_n19_r06_all_flat = rm_n19_r06_all.flatten()
    vgac_r06_all_flat = vgac_r06_all.flatten()
    viirs_r06_all_flat = viirs_r06_all.flatten()
    
#     fig = plt.figure()
#     fig.suptitle('r06%s' %title_end)
#     ax = fig.add_subplot(1,2,1)
#     ax.scatter(rm_n19_r06_all_flat, vgac_r06_all_flat)
#     ax.set_title('Scatter')
# #     ax.set_title('r06%s' %title_end)
#     ax.set_xlabel('NOAA 19, max = %d' %int(rm_n19_r06_all_flat.max()))
#     ax.set_ylabel('VGAC, max = %d' %int(vgac_r06_all_flat.max()))
#     ax.set_xlim([0, r_max_axis])
#     ax.set_ylim([0, r_max_axis])
#     ax.set_xticks(r_plot_ticks_int)
#     ax.set_yticks(r_plot_ticks_int)
#     
#     xn = np.linspace(0, r_max_axis, 100)
#     k2, m = np.ma.polyfit(rm_n19_r06_all_flat, vgac_r06_all_flat, 1)
#     ny_y2 = k2*xn + m #rm_n19_r06_all_flat + m
# #     ax.plot(rm_n19_r06_all_flat, ny_y2, 'g*', label = '%.2G * x %+.2f' %(k2, m))
#     ax.plot(xn, ny_y2, 'g--', label = '%.4G*x%+.2f' %(k2, m))
# #     inds = np.where(~(rm_n19_r06_all_flat.mask | vgac_r06_all_flat.mask))
# #     rp = scipy.stats.pearsonr(rm_n19_r06_all_flat[inds], vgac_r06_all_flat[inds])
# #     rs = scipy.stats.spearmanr(rm_n19_r06_all_flat[inds], vgac_r06_all_flat[inds])
#     rr = np.ma.corrcoef(rm_n19_r06_all_flat, vgac_r06_all_flat)
#     
# #     a_r06, b_r06 = np.polyfit(rm_n19_r06.flatten(), vgac_r06.flatten(), 1)
# #     ax.plot(xn, (a_r06 + (b_r06 * xn)), 'r')
#     ax.plot(xn, xn, 'r', label = '1 * x')
#     ax.text(1,0,'ccof = %.4f' %rr[0,1], horizontalalignment='right', transform=ax.transAxes)#, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
#     ax.legend(loc=2)
#     ax.set_aspect(1)
# 
#     ax = fig.add_subplot(1,2,2)
# #     ax.hist2d(rm_n19_r06_all_flat, vgac_r06_all_flat)#, bins=histd2_bins, range=[0, r_max_axis] )
# #     ax = fig.add_subplot(1,3,3)
#     H, xedges, yedges = np.histogram2d(rm_n19_r06_all_flat, vgac_r06_all_flat, bins=histd2_bins, range=[[0, r_max_axis], [0, r_max_axis]])
#     cmap = 'plasma'
#     im = ax.imshow(H.T, origin='lower', vmin=0, vmax=vmax, cmap=cmap)
#     #: -1 since the lat one is included in contrast to range. -0.5 to get everything inside the plot
#     xn = np.linspace(0, histd2_bins-1.5, 100)
# #     xn = range(histd2_bins)
#     ny_y2 = k2*xn + m 
#     ax.plot(xn, ny_y2, 'g--')#, label = '%.4G*x%+.2f' %(k2, m))
#     ax.plot(xn, xn, 'r')
# #     ax.set_xticks([0, 25, 50, 75, 99])
# #     ax.set_xticklabels(['0', '45', '90', '135', '180'])
#     ax.set_title('2D Histogram')
#     ax.set_xticks(r_plot_ticks_loc)
#     ax.set_xticklabels(r_plot_ticks_str)
#     ax.set_yticks(r_plot_ticks_loc)
#     ax.set_yticklabels(r_plot_ticks_str)
#     fig.subplots_adjust(right=0.89)
#     pos2 = ax.get_position()
#     cbar_ax = fig.add_axes([0.90, pos2.y0, 0.01, pos2.y1 - pos2.y0])
#     cbar = fig.colorbar(im, cax=cbar_ax)
#     # fig.tight_layout(rect=[0, 0.03, 0.97, 0.97])
#     # fig.tight_layout(rect=[0, 0.03, 1, 0.95])
#     # plt.tight_layout()
# #     ax.set_xticks(xedges)
# #     fig.tight_layout()
#     figname = '%s/%s_n19_vgac_r06_satz-%d_sunz-%d_td-%dmin' %(plotDir, plotType, accept_satz_max, accept_sunz_max, accept_time_diff)
#     fig.savefig(figname + '.png')
#     fig.show()
    
    
    #: VGAC
    figname = '%s/%s_n19_vgac_r06_satz-%d_sunz-%d_td-%dmin' %(plotDir, plotType, accept_satz_max, accept_sunz_max, accept_time_diff)
    plotScatterHisto(rm_n19_r06_all_flat, vgac_r06_all_flat, 'r06%s' %title_end, 'VGAC', [0, r_max_axis], \
                     r_plot_ticks_int, r_plot_ticks_str, r_plot_ticks_loc, histd2_bins, vmax, figname)
    #:VIIRS
    figname = '%s/%s_n19_viirs_r06_satz-%d_sunz-%d_td-%dmin' %(plotDir, plotType, accept_satz_max, accept_sunz_max, accept_time_diff)
    plotScatterHisto(rm_n19_r06_all_flat, viirs_r06_all_flat, 'r06%s' %title_end, 'VIIRS', [0, r_max_axis], \
                     r_plot_ticks_int, r_plot_ticks_str, r_plot_ticks_loc, histd2_bins, vmax, figname)
    
    #: Channel r09
    rm_n19_r09_all_flat = rm_n19_r09_all.flatten()
    vgac_r09_all_flat = vgac_r09_all.flatten()
    viirs_r09_all_flat = viirs_r09_all.flatten()
    #: VGAC
    figname = '%s/%s_n19_vgac_r09_satz-%d_sunz-%d_td-%dmin' %(plotDir, plotType, accept_satz_max, accept_sunz_max, accept_time_diff)
    plotScatterHisto(rm_n19_r09_all_flat, vgac_r09_all_flat, 'r09%s' %title_end, 'VGAC', [0, r_max_axis], \
                     r_plot_ticks_int, r_plot_ticks_str, r_plot_ticks_loc, histd2_bins, vmax, figname)

    #: VIIRS
    figname = '%s/%s_n19_viirs_r09_satz-%d_sunz-%d_td-%dmin' %(plotDir, plotType, accept_satz_max, accept_sunz_max, accept_time_diff)
    plotScatterHisto(rm_n19_r09_all_flat, viirs_r09_all_flat, 'r09%s' %title_end, 'VIIRS', [0, r_max_axis], \
                     r_plot_ticks_int, r_plot_ticks_str, r_plot_ticks_loc, histd2_bins, vmax, figname)
    
    #: Channel tb11
    rm_n19_tb11_all_flat = rm_n19_tb11_all.flatten()
    vgac_tb11_all_flat = vgac_tb11_all.flatten()
    viirs_tb11_all_flat = viirs_tb11_all.flatten()
#     if plotType == 'density':
#         print('Calculate Density')
#         tic = time.time()
#     
#         xy = np.ma.vstack([rm_n19_tb11_all.flatten(), vgac_tb11_all.flatten()])
#         kde = gaussian_kde(xy.data[~xy.mask])
#         z = kde(xy.data[~xy.mask])
#         
#         idx = z.argsort()
#     
#         rm_n19_tb11_all_flat, vgac_tb11_all_flat, z = rm_n19_tb11_all_flat[idx], vgac_tb11_all_flat[idx], z[idx]
#         cc=z
#         ss=50
#         print(time.time() - tic)
#         pdb.set_trace()
#     else:
#         cc=None
#         ss=None
    #: VGAC
    figname = '%s/%s_n19_vgac_tb11_satz-%d_sunz-%d_td-%dmin' %(plotDir, plotType, accept_satz_max, accept_sunz_max, accept_time_diff) 
    plotScatterHisto(rm_n19_tb11_all_flat, vgac_tb11_all_flat, 'tb11%s' %title_end, 'VGAC', [tb_min_axis, tb_max_axis], \
                     tb_plot_ticks_int, tb_plot_ticks_str, tb_plot_ticks_loc, histd2_bins, vmax, figname)
    #: VIIRS
    figname = '%s/%s_n19_viirs_tb11_satz-%d_sunz-%d_td-%dmin' %(plotDir, plotType, accept_satz_max, accept_sunz_max, accept_time_diff)
    plotScatterHisto(rm_n19_tb11_all_flat, viirs_tb11_all_flat, 'tb11%s' %title_end, 'VIIRS', [tb_min_axis, tb_max_axis], \
                     tb_plot_ticks_int, tb_plot_ticks_str, tb_plot_ticks_loc, histd2_bins, vmax, figname)
    
    #: Channel tb12
    rm_n19_tb12_all_flat = rm_n19_tb12_all.flatten()
    vgac_tb12_all_flat = vgac_tb12_all.flatten()
    viirs_tb12_all_flat = viirs_tb12_all.flatten()
    #: VGAC
    
    figname = '%s/%s_n19_vgac_tb12_satz-%d_sunz-%d_td-%dmin' %(plotDir, plotType, accept_satz_max, accept_sunz_max, accept_time_diff)
    plotScatterHisto(rm_n19_tb12_all_flat, vgac_tb12_all_flat, 'tb12%s' %title_end, 'VGAC', [tb_min_axis, tb_max_axis], \
                     tb_plot_ticks_int, tb_plot_ticks_str, tb_plot_ticks_loc, histd2_bins, vmax, figname)
    #: VIIRS
    figname = '%s/%s_n19_viirs_tb12_satz-%d_sunz-%d_td-%dmin' %(plotDir, plotType, accept_satz_max, accept_sunz_max, accept_time_diff)
    plotScatterHisto(rm_n19_tb12_all_flat, viirs_tb12_all_flat, 'tb12%s' %title_end, 'VIIRS', [tb_min_axis, tb_max_axis], \
                     tb_plot_ticks_int, tb_plot_ticks_str, tb_plot_ticks_loc, histd2_bins, vmax, figname)

    #: Channel tb37
    rm_n19_tb37_all_flat = rm_n19_tb37_all.flatten()
    vgac_tb37_all_flat = vgac_tb37_all.flatten()
    viirs_tb37_all_flat = viirs_tb37_all.flatten()
    #: VGAC
    figname = '%s/%s_n19_vgac_tb37_satz-%d_sunz-%d_td-%dmin' %(plotDir, plotType, accept_satz_max, accept_sunz_max, accept_time_diff)
    plotScatterHisto(rm_n19_tb37_all_flat, vgac_tb37_all_flat, 'tb37%s' %title_end, 'VGAC', [tb_min_axis, tb_max_axis], \
                     tb_plot_ticks_int, tb_plot_ticks_str, tb_plot_ticks_loc, histd2_bins, vmax, figname)
    
    #: VIIRS
    figname = '%s/%s_n19_viirs_tb37_satz-%d_sunz-%d_td-%dmin' %(plotDir, plotType, accept_satz_max, accept_sunz_max, accept_time_diff) 
    plotScatterHisto(rm_n19_tb37_all_flat, viirs_tb37_all_flat, 'tb37%s' %title_end, 'VIIRS', [tb_min_axis, tb_max_axis], \
                     tb_plot_ticks_int, tb_plot_ticks_str, tb_plot_ticks_loc, histd2_bins, vmax, figname)
    
    
    
    pdb.set_trace()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    