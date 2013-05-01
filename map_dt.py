#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  map_dt.py
#   Purpose:   maps dt calculated with FFM
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
import os
import sys

# ------------------- INPUT -----------------------------
processed_events_add = '/import/neptun-radler/hosseini-downloads/KASRA/FFM'
band = 'band01'
band = 'BB'
xcorr_limit = 0.9
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
passed_staev = []
evlats = []
evlons = []
print '\nERRORS:'
for i in range(len(proc_ev_ls)):
    evnt = proc_ev_ls[i]
    all_dt_event = np.array([])
    passed_staev_tmp = []
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
                passed_staev_tmp.append([lat, lon, xcorr, float(evlat), float(evlon), i])
                all_dt_event = np.append(all_dt_event, dt)
        all_dt_median = all_dt_event - np.median(all_dt_event)
        for k in range(len(all_dt_median)):
            passed_staev_tmp[k].insert(2, all_dt_median[k])
            passed_staev.append(passed_staev_tmp[k])
        evlats.append(float(evlat))
        evlons.append(float(evlon))
    except Exception, e:
        print e
        failed += 1

print '\n------------------------'
xcorr_sum = 0
dt_sum = 0
for i in range(len(passed_staev)):
    xcorr_sum += passed_staev[i][3]
    dt_sum += passed_staev[i][2]

print '%s events failed!' %(failed)
print '%s station-event pairs found...' %(len(passed_staev))
print 'mean xcorr: %s' %(xcorr_sum/len(passed_staev))
print 'mean dt: %s' %(dt_sum/len(passed_staev))
print '------------------------'

print '\nPlotting...'

lon_inc = [-180, -90, 0, 90, 90, 90]
lat_inc = [0, 0, 0, 0, 90, -90]
for inc in range(len(lon_inc)):
    print '\nPlot for lon: %s and lat: %s view' %(lon_inc[inc], lat_inc[inc])
    plt.subplot(3,2,inc+1)
    #m = Basemap(projection='aeqd', lon_0=-100, lat_0=0, resolution='c')
    #m = Basemap(projection='cyl', lon_0=0, lat_0=0, resolution='c')
    if not lat_inc[inc] in [90, -90]:
        m = Basemap(projection='moll', lon_0=lon_inc[inc], lat_0=lat_inc[inc], resolution='c')
    elif lat_inc[inc] == 90:
        m = Basemap(projection='npaeqd',boundinglat=30,lon_0=270,resolution='l')
    else:
        m = Basemap(projection='spaeqd',boundinglat=-30,lon_0=270,resolution='l')
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
                    passed_staev[i][1], passed_staev[i][0], int(divisions))
        gc_xmid_points.append(gc_arc[0][int(divisions/2)])
        gc_ymid_points.append(gc_arc[1][int(divisions/2)])
        gc_dtmid_points.append(passed_staev[i][2])
        #m.scatter(gc_arc[0][50], gc_arc[1][50], color=passed_staev[i][-1])
        #m.scatter(gc_arc[0][50], gc_arc[1][50], color=passed_staev[i][2], vmin=-6, vmax=6.0)

    map_evlats = []
    map_evlons = []
    print '\nError:' 
    for i in range(len(proc_ev_ls)):
        try:
            evx, evy = m(evlons[i], evlats[i])
            map_evlons.append(evx)
            map_evlats.append(evy)
        except Exception, e:
            print e

    # XXX just to get a better figure!
    gc_xmid_points_sub = []
    gc_ymid_points_sub = []
    gc_dt_median_sub = []
    for i in range(len(gc_dtmid_points)):
        if -4<=gc_dtmid_points[i]<=4:
            gc_dt_median_sub.append(gc_dtmid_points[i])
            gc_xmid_points_sub.append(gc_xmid_points[i])
            gc_ymid_points_sub.append(gc_ymid_points[i])

    m.scatter(gc_xmid_points_sub, gc_ymid_points_sub, c=gc_dt_median_sub,
                vmin=min(gc_dt_median_sub), vmax=max(gc_dt_median_sub), 
                alpha=0.6, edgecolor='none', zorder=10)

    m.colorbar()
    m.scatter(map_evlons, map_evlats, s=120, c='r', edgecolor='none', zorder=10, marker='*')
    plt.title('Dominant Period: %s\n#event-station pairs: %s' %(band_p, len(passed_staev)))
plt.show()
print '\n'

# XXX Still not sure about the colors!!!!!!!
#for i in range(len(gc_dt_median)):
#    gc_dt_median[i] = 1*gc_dt_median[i]
#m.scatter(gc_xmid_points, gc_ymid_points, c=gc_dt_median, vmin=-3, vmax=3, 
#            alpha=0.6, edgecolor='none', zorder=10)

# XXX just to get a better figure!
#gc_xmid_points_new = []
#gc_ymid_points_new = []
#gc_dt_median_new = []
#for i in range(len(gc_dt_median)):
#    if -4<=gc_dt_median[i]<=4:
#        gc_dt_median_new.append(gc_dt_median[i])
#        gc_xmid_points_new.append(gc_xmid_points[i])
#        gc_ymid_points_new.append(gc_ymid_points[i])
#
#m.scatter(gc_xmid_points_new, gc_ymid_points_new, c=gc_dt_median_new, vmin=min(gc_dt_median_new), vmax=max(gc_dt_median_new), 
#            alpha=0.6, edgecolor='none', zorder=10)

#print 'Median: %s' %(np.median(gc_dtmid_points))
#gc_dt_median = gc_dtmid_points - np.median(gc_dtmid_points)
#m.scatter(gc_xmid_points, gc_ymid_points, c=gc_dt_median, vmin=min(gc_dt_median), vmax=max(gc_dt_median), 
#            alpha=0.6, edgecolor='none', zorder=10)

#m.scatter(gc_xmid_points, gc_ymid_points, c=gc_dtmid_points, vmin=min(gc_dtmid_points), vmax=max(gc_dtmid_points), alpha=0.6, edgecolor='none', zorder=10)
