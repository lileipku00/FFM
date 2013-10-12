#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  plot_dT_dA.py
#   Purpose:   plot dT and dA againts dominant period for FFM
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import glob
import os
import shutil
import sys

import util_ffproc as uf 

# ------------------- INPUT -----------------------------
xcorr_limit = 0.8
all_processed = True
update_all = True
plot_scatter = False
plot_mean = True
plot_sta = True
all_stations = True
just_selected_events = False
list_events = [\
#'0007.2009.089.a',  
#'0013.2009.136.a',
#'0022.2009.321.a',
#'0049.2009.184.a',
#'0054.2009.267.a',
#'0070.2009.123.a',
#'0097.2009.255.a',
#'0104.2009.130.a',
#'0109.2009.085.a',
#'0115.2009.033.a',
#'0118.2009.273.a',
#'0122.2009.107.a',
#'0153.2009.059.a',
#'0200.2009.091.a',
#'0211.2009.224.a',
#'0221.2009.253.a',
#'0222.2009.015.a',
#'0224.2009.156.a',
#'0228.2009.157.a',
#'0230.2009.221.a',
#'0234.2009.085.a',
#'0235.2009.246.a',
#'0238.2009.271.a',
#'0250.2009.261.a',
#'0263.2009.053.a',
#'0265.2009.251.a',
#'0267.2009.261.a',
#'0274.2009.273.a',
#'0188.2009.015.a',
#'0181.2009.026.a',
#'0167.2009.266.a',
#'0169.2009.272.a',
'0161.2009.196.a',
'0161.2009.217.a',
'0171.2009.230.a',
'0173.2009.065.a',
'0173.2009.211.a',
'0174.2009.078.a',
'0177.2009.144.a',
'0178.2009.018.a',
#'0190.2009.174.a',
'0192.2009.022.a',
'0196.2009.003.a',
'0284.2009.261.a',
'0312.2009.264.a',
'0325.2009.240.a',
'0370.2009.182.c',
'0381.2009.096.a',
'0403.2009.157.a',
'0412.2009.167.a',
'0429.2009.085.a',
'0429.2009.132.a',
'0431.2009.155.a',
'0641.2009.065.a',
'0681.2009.188.a',
'0684.2009.001.a',
'0685.2009.260.a',
'0691.2009.075.a',
'0696.2009.074.a',
'0702.2009.094.a',
'0703.2009.222.a',
'0718.2009.302.a',
'0731.2009.148.a',
'0756.2009.189.a']

remote_dir = '/import/neptun-helles/hosseini/FFM'
# -------------------------------------------------------

intro = 20*'-'
intro += '\nUsage:'
intro += '\npython plot_dT_dA.py bands[1-3] address-for-FFM-event.\n'
intro += 20*'-'

bands = sys.argv[1]
bands = range(int(bands[0]), int(bands[-1])+1)
band_period = {'1': 30.0,'2': 21.2,'3': 15.0,'4': 10.6,'5': 7.5,
                '6': 5.3,'7': 3.7,'8': 2.7}
if not all_processed:
    evadd = sys.argv[2]
    evname = evadd.split('/')[-1]
    if not evname: evname = evadd.split('/')[-2]
    all_staev = uf.reader(evadd, bands, band_period, all_stations=all_stations)
    passed_staev = uf.filters(all_staev, bands, xcorr_limit=xcorr_limit,
                                                all_stations=all_stations)
    if plot_scatter: uf.ffpscatter(passed_staev)
    if plot_mean: 
        per, dt_mean, a_mean, flag = uf.stamean(passed_staev)
        uf.ffplotmean(per, dt_mean)
else:
    if not just_selected_events:
        proc_ev_ls = glob.glob(os.path.join(remote_dir, '*.*.*.*'))
    else:
        proc_ev_ls = []
        for ev in list_events:
            proc_ev_ls.append(os.path.join(remote_dir, ev))
    print '%s processed events found!' %(len(proc_ev_ls))
    if update_all:
        shutil.rmtree(os.path.join('.', 'statistics'))
        os.mkdir(os.path.join('.', 'statistics'))
    all_passed_staev = []
    for i in range(len(proc_ev_ls)):
        all_staev = uf.reader(proc_ev_ls[i], bands, band_period, all_stations=all_stations)
        if all_staev == []: continue
        passed_staev = uf.filters(all_staev, bands, xcorr_limit=xcorr_limit, 
                                                all_stations=all_stations)
        if passed_staev[0] == []: continue
        all_passed_staev.append(passed_staev)
        if update_all: uf.writer(passed_staev, bands)
    if plot_scatter: uf.ffpscatter(all_passed_staev, all_events=True)
    if plot_mean:
        all_dt_mean = []
        all_a_mean = []
        all_lat = []
        all_lon = []
        for i in range(len(all_passed_staev)):
            per, dt_mean, a_mean, flag = uf.stamean(all_passed_staev[i])
            if flag:
                all_dt_mean.append([dt_mean, len(all_passed_staev[i][0])])
                all_a_mean.append([a_mean, len(all_passed_staev[i][0])])
                for sta in all_passed_staev[i][0]:
                    all_lat.append(sta[0])
                    all_lon.append(sta[1])
        #uf.allffplotmean(per, all_dt_mean)
        uf.meanall_ffplot(per, all_dt_mean, all_a_mean)

    if plot_sta:
        import matplotlib.pyplot as plt
        from mpl_toolkits.basemap import Basemap
        import numpy as np
        m = Basemap(projection='cyl', lon_0=0.0, lat_0=0.0, resolution='c')
        m.drawcoastlines()
        m.fillcontinents()
        m.drawparallels(np.arange(-90., 120., 30.))
        m.drawmeridians(np.arange(0., 420., 60.))
        m.drawmapboundary()
        x, y = m(all_lon, all_lat)
        m.scatter(x, y, c='red', edgecolor='none', zorder=20, marker='v', s=40)
        plt.show()
