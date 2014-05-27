#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  find_events_azimuth.py
#   Purpose:   Adapted from find_events to find interesting events!
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python and Obspy modules will be imported in this part.
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from obspy.core.util import gps2DistAzimuth
from obspy.imaging.beachball import Beach
import os

# ------------------- INPUT -----------------------------
# min and max desired magnitude:
min_mag = 6.0
max_mag = 10.0

# rectangular window to select events:
min_lat = -90.
max_lat = 90.
min_lon = -180.0
max_lon = 180.0
# Around Hawaii:
#min_lat = 13.7
#max_lat = 28.
#min_lon = -167.2
#max_lon = -148.5

# First we define a center point (center_lat and center_lon)
# Pick an azimuth which is interesting for you
# Program searches for those events that they are within a range (azimuth_error) from that specific point
center_lat = 40.
center_lon = -115.
req_azimuth = 312.37
azimuth_error = 3.0

# min and maximum distance to the center point
min_dist = 32
max_dist = 90

min_date = 2008
max_date = 2014

plot_ev = True
# Center point to plot the "nsper" projection
center_lat_plot = -0.72+40
center_lon_plot = 99.867+60

# First lat and lon to plot the great circle distance
# This is usually the coordinates of the event that has been manually picked
event_lat_gc = -0.72
event_lon_gc = 99.867
# -------------------------------------------------------

fio_pdata_events = open(os.path.join(os.path.curdir, 'results', 'pdata_events.txt'), 'r')
pdata_events = fio_pdata_events.readlines()[1:]
print '%s total number of events' % len(pdata_events)
fio_selected_events = open(os.path.join('.', 'results', 'selected_events.txt'), 'w')

m2degree = 360./(2.*np.pi*6371000.)
enum = 0
for i in range(len(pdata_events)):

    if not min_mag <= float(pdata_events[i].split(',')[4]) <= max_mag: continue

    ev_lat = float(pdata_events[i].split(',')[1])
    ev_lon = float(pdata_events[i].split(',')[2])

    if not min_lat <= ev_lat <= max_lat: continue
    if not min_lon <= ev_lon <= max_lon: continue

    cent_ev_gd = gps2DistAzimuth(center_lat, center_lon, ev_lat, ev_lon)

    if not abs(cent_ev_gd[1] - req_azimuth) <= azimuth_error: continue
    if not min_date <= int(pdata_events[i].split(',')[0].split('.')[1]) <= max_date: continue

    dist = cent_ev_gd[0]*m2degree
    if not min_dist <= dist <= max_dist: continue

    fio_selected_events.writelines(pdata_events[i])
    print '-------------'
    print pdata_events[i].split(',')[0]
    print 'info:\nDistance: %s' % dist
    print pdata_events[i]
    enum += 1

print '\n===================='
print 'PASSED events: %s' % enum
print '===================='
fio_selected_events.close()

plt.ion()
if plot_ev:
    print 'Plotting the events...'
    #m = Basemap(projection='robin', lon_0=99.867, lat_0=0.0, resolution='c')
    h = 300000000
    m = Basemap(projection='nsper', lon_0=center_lon_plot, lat_0=center_lat_plot, satellite_height=h*1000.,
                resolution='l')
    m.drawcoastlines()
    m.fillcontinents()
    m.drawparallels(np.arange(-90., 120., 30.))
    m.drawmeridians(np.arange(0., 420., 60.))
    m.drawmapboundary()
    #m.drawgreatcircle(event_lon_gc, event_lat_gc, center_lon, center_lat, linewidth=3)

    fio_selected_events = open(os.path.join(os.path.curdir, 'results', 'selected_events.txt'), 'r')
    selected_events = fio_selected_events.readlines()
    for i in range(len(selected_events)):
        try:
            evnt = selected_events[i].split(',')
            x, y = m(float(evnt[2]), float(evnt[1]))
            #m.scatter(x, y, c='b', s=400, marker='v', zorder=200)
            m.drawgreatcircle(float(evnt[2]), float(evnt[1]), center_lon, center_lat, linewidth=3)
            focmecs = [float(evnt[5]), float(evnt[6]), float(evnt[7]), float(evnt[8]), float(evnt[9]), float(evnt[10])]
            ax = plt.gca()
            b = Beach(focmecs, xy=(x, y), width=7e5, linewidth=1, alpha=0.85, size=400)
            b.set_zorder(10)
            ax.add_collection(b)
        except Exception, e:
            print 'ERROR in %s:\n%s' % (evnt[0], e)


    ##  # ------------------------------------------------
    ##  # MANUALLY CREATE A FIGURE:
    ##  center_lat = 40.
    ##  center_lon = -115.

    ##  # EVENT-1:
    ##  # 0274.2009.273.a,-0.7200000,99.86700,82.00000,7.5,0.163E+21,-0.148E+20,-0.148E+21,0.349E+20,-0.397E+20,-0.157E+21
    ##  event_lat_gc = -0.7200000
    ##  event_lon_gc = 99.86700

    ##  x, y = m(event_lon_gc, event_lat_gc)
    ##  focmecs = [0.163E+21, -0.148E+20, -0.148E+21, 0.349E+20, -0.397E+20, -0.157E+21]
    ##  ax = plt.gca()
    ##  b = Beach(focmecs, xy=(x, y), width=7e5, linewidth=1, alpha=0.85, size=400)
    ##  b.set_zorder(10)
    ##  ax.add_collection(b)

    ##  # EVENT-2
    ##  # '0221.2009.111.a', '50.83300', '155.0090', '152.0000', '6.2',
    ##  x, y = m(155.0090, 50.83300)
    ##  focmecs = [-0.229E+18, -0.368E+18, 0.596E+18, -0.130E+19, -0.201E+19, 0.194E+18]
    ##  ax = plt.gca()
    ##  b = Beach(focmecs, xy=(x, y), width=7e5, linewidth=1, alpha=0.85, size=400)
    ##  b.set_zorder(10)
    ##  ax.add_collection(b)

    ##  m.drawgreatcircle(155.0090, 50.83300, center_lon, center_lat, linewidth=3, c='r')
    ##  m.drawgreatcircle(event_lon_gc, event_lat_gc, center_lon, center_lat, linewidth=3, c='b')

    ##  ### # EVENT-3
    ##  ### # 0742.2008.348.a,-48.90000,123.3000,2.000000,5.9,-0.115E+18,0.687E+18,-0.572E+18,-0.744E+18,0.615E+17,-0.206E+18,
    ##  ### event_lat_gc = -48.90
    ##  ### event_lon_gc = 123.30

    ##  ### x, y = m(event_lon_gc, event_lat_gc)
    ##  ### focmecs = [-0.115E+18, 0.687E+18, -0.572E+18, -0.744E+18, 0.615E+17, -0.206E+18]
    ##  ### ax = plt.gca()
    ##  ### b = Beach(focmecs, xy=(x, y), width=7e5, linewidth=1, alpha=0.85, size=400)
    ##  ### b.set_zorder(10)
    ##  ### ax.add_collection(b)

    ##  ### # EVENT-4
    ##  ### # 0613.2007.226.a,19.34900,-155.0730,8.000000,5.6,0.745E+17,-0.138E+17,-0.608E+17,0.969E+17,0.822E+17,-0.592E+17,
    ##  ### x, y = m(-155.0730, 19.34900)
    ##  ### focmecs = [0.745E+17, -0.138E+17, -0.608E+17, 0.969E+17, 0.822E+17, -0.592E+17]
    ##  ### ax = plt.gca()
    ##  ### b = Beach(focmecs, xy=(x, y), width=7e5, linewidth=1, alpha=0.85, size=400)
    ##  ### b.set_zorder(10)
    ##  ### ax.add_collection(b)

    ##  ### m.drawgreatcircle(-155.0730, 19.34900, center_lon, center_lat, linewidth=3, c='r')
    ##  ### m.drawgreatcircle(event_lon_gc, event_lat_gc, center_lon, center_lat, linewidth=3, c='b')

    ##  x, y = m(center_lon, center_lat)
    ##  m.scatter(x, y, c='r', s=600, linewidth=5, marker='x', zorder=200)

    ##  plt.show()
