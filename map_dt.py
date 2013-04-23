#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  map_dt.py
#   Purpose:   plot dt maps based on the Finite Frequency results 
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python and Obspy modules will be imported in this part.
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import os
import sys

# ------------------- INPUT -----------------------------
processed_events_add = '/import/neptun-radler/hosseini-downloads/KASRA/FFM'
band = 'band01'
#band = 'BB'
xcorr_limit = 0.9
proc_ev_ls = glob.glob(os.path.join(processed_events_add, '*.*.*.*'))
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
else:
    sys.exit('Wrong band!')

print '%s processed events found!' %(len(proc_ev_ls))
failed = 0
passed_staev = []
evlats = []
evlons = []
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
        f_dt = fio_dt.readlines()
        for j in range(2, len(f_dt)):
            info_dt = f_dt[j].split()
            xcorr = float(info_dt[2])
            dt = float(info_dt[5])
            lat = float(info_dt[6])
            lon = float(info_dt[7])
            if xcorr >= xcorr_limit:
                passed_staev.append([lat, lon, dt, xcorr, float(evlat), float(evlon)])
        evlats.append(float(evlat))
        evlons.append(float(evlon))
    except Exception, e:
        print e
        failed += 1

xcorr_sum = 0
dt_sum = 0
for i in range(len(passed_staev)):
    xcorr_sum += passed_staev[i][3]
    dt_sum += passed_staev[i][2]

print '%s events failed!' %(failed)
print '%s station-event found...' %(len(passed_staev))
print 'mean xcorr: %s' %(xcorr_sum/len(passed_staev))
print 'mean dt: %s' %(dt_sum/len(passed_staev))

print 'Plotting...'
#m = Basemap(projection='aeqd', lon_0=-100, lat_0=40, resolution='c')
m = Basemap(projection='cyl', lon_0=0, lat_0=0, resolution='c')
m.drawcoastlines()
#m.fillcontinents()
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
m.drawmapboundary()
gc_xmid_points = []
gc_ymid_points = []
gc_dtmid_points = []
for i in range(len(passed_staev)):
    sys.stdout.write('\r')
    sys.stdout.write("[%-100s] %d%%" % ('='*int(100.*(i+1)/len(passed_staev)),
                                            100.*(i+1)/len(passed_staev)))
    sys.stdout.flush()
    #m.drawgreatcircle(passed_staev[i][5],passed_staev[i][4],
    #            passed_staev[i][1],passed_staev[i][0], alpha = 0.1, 
    #            color=passed_staev[i][-1])
    gc_arc = m.gcpoints(passed_staev[i][5], passed_staev[i][4],
                passed_staev[i][1], passed_staev[i][0], 100)
    gc_xmid_points.append(gc_arc[0][50])
    gc_ymid_points.append(gc_arc[1][50])
    gc_dtmid_points.append(passed_staev[i][2])
    #m.scatter(gc_arc[0][50], gc_arc[1][50], color=passed_staev[i][-1])
    #m.scatter(gc_arc[0][50], gc_arc[1][50], color=passed_staev[i][2], vmin=-6, vmax=6.0)

map_evlats = []
map_evlons = []
for i in range(len(proc_ev_ls)):
    try:
        evx, evy = m(evlons[i], evlats[i])
        map_evlons.append(evx)
        map_evlats.append(evy)
    except Exception, e:
        print e
gc_dt_median = gc_dtmid_points - np.median(gc_dtmid_points)
# XXX Still not sure about the colors!!!!!!!
for i in range(len(gc_dt_median)):
    gc_dt_median[i] = -1*gc_dt_median[i]
m.scatter(gc_xmid_points, gc_ymid_points, c=gc_dt_median, vmin=-3, vmax=3, 
            alpha=0.6, edgecolor='none', zorder=10)
m.colorbar()
m.scatter(map_evlons, map_evlats, s=120, c='r', edgecolor='none', zorder=10, marker='*')
plt.title('Dominant Period: %s\n#event-station pairs: %s' %(band_p, len(passed_staev)))
plt.show()
print '\n Ciao!'
