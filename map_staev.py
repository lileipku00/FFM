#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  map_staev.py
#   Purpose:   maps stations and events used in measurements
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from obspy.imaging.beachball import Beach
import os
import sys

# ------------------- INPUT -----------------------------
processed_events_add = '/import/neptun-helles/hosseini/FFM'
band = 'band01'
#band = 'BB'
xcorr_limit = 0.85
# Number of divisions on gcarc for plotting reasons (find the middle point)
divisions = 3
# -------------------------------------------------------

band_period = {
    'band01': '30.0sec',
    'band02': '21.2sec',
    'band03': '15.0sec',
    'band04': '10.6sec',
    'band05': '7.5sec',
    'band06': '5.3sec',
    'band07': '3.7sec',
    'band08': '2.7sec'}

if band in band_period:
    band_p = band_period[band]
elif band == 'BB':
    band_p = 'BB'
else:
    sys.exit('Wrong band is specified! (check the INPUT)')

proc_ev_ls = glob.glob(os.path.join(processed_events_add, '*.*.*.*'))
print '%s processed events found!' %(len(proc_ev_ls))
failed = 0
events_info = []
stations_info = []
print '\nERRORS:'
for i in range(len(proc_ev_ls)):
    evnt = proc_ev_ls[i]
    try:
        fio_dt = open(os.path.join(evnt, 'outfiles', 'ffproc.ampstt.' + band), 'r')
        fio_source = open(os.path.join(evnt, 'outfiles', 'ampinv.source'), 'r')
        f_source = fio_source.readlines()
        ev_year, ev_julianday, ev_hr, ev_min, ev_sec, ev_msec = f_source[1].split()
        evlat, evlon, catalog_depth, inverted_depth = f_source[3].split()
        try:
            mrr, mtt, mpp, mrt, mrp, mtp = f_source[13].split()
        except Exception, e:
            mrr, mtt, mpp, mrt, mrp, mtp = f_source[7].split()
        events_info.append([float(evlat), float(evlon), mrr, mtt, mpp, mrt, mrp, mtp])
        f_dt = fio_dt.readlines()
        for j in range(2, len(f_dt)):
            info_dt = f_dt[j].split()
            xcorr = float(info_dt[2])
            lat = float(info_dt[6])
            lon = float(info_dt[7])
            if xcorr >= xcorr_limit:
                stations_info.append([lat, lon])
    except Exception, e:
        print e
        failed += 1

print '%s events failed!' %(failed)
print '%s station-event pairs found...' %(len(stations_info))
print '------------------------'

stations_info_trim = []
for item in stations_info:
    if item not in stations_info_trim:
        stations_info_trim.append(item)

print '\nPlotting...'

m = Basemap(projection='cyl', lon_0=0.0, lat_0=0.0, resolution='c')
m.drawcoastlines()
m.fillcontinents()
m.drawparallels(np.arange(-90., 120., 30.))
m.drawmeridians(np.arange(0., 420., 60.))
m.drawmapboundary()
for i in range(len(events_info)):
    print i,
    sys.stdout.flush()
    try:
        x, y = m(float(events_info[i][1]), float(events_info[i][0]))
        focmecs = [float(events_info[i][2]),float(events_info[i][3]),float(events_info[i][4]),
                    float(events_info[i][5]),float(events_info[i][6]),float(events_info[i][7])]
        ax = plt.gca()
        b = Beach(focmecs, xy=(x, y), width=3, linewidth=1, alpha=0.85)
        b.set_zorder(10)
        ax.add_collection(b)
    except Exception, e:
        print evnt[0], e

for i in range(len(stations_info_trim)):
    print i,
    sys.stdout.flush()
    try:
        x, y = m(float(stations_info[i][1]), float(stations_info[i][0]))
        m.scatter(x, y, c='red', edgecolor='none', zorder=20, marker='v', s=40)
    except Exception, e:
        print evnt[0], e

plt.show()
