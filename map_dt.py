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
xcorr_limit = 0.95
proc_ev_ls = glob.glob(os.path.join(processed_events_add, '*.*.*.*'))
# -------------------------------------------------------

print '%s processed events found!' %(len(proc_ev_ls))
failed = 0
passed_staev = []
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
                # XXX not sure about the colors!!!
                if dt >= 0: color = 'r'
                else: color = 'b'
                passed_staev.append([lat, lon, dt, xcorr, float(evlat), float(evlon), color])
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
m = Basemap(projection='aeqd', lon_0=-100, lat_0=40, resolution='c')
m.drawcoastlines()
m.fillcontinents()
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
m.drawmapboundary()

for i in range(len(passed_staev)):
    sys.stdout.write('\r')
    sys.stdout.write("[%-100s] %d%%" % ('='*int(100.*(i+1)/len(passed_staev)),
                                            100.*(i+1)/len(passed_staev)))
    sys.stdout.flush()
    m.drawgreatcircle(passed_staev[i][5],passed_staev[i][4],
                passed_staev[i][1],passed_staev[i][0], alpha = 0.1, 
                color=passed_staev[i][-1])

plt.show()
