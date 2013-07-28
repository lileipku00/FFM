#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  find_events.py
#   Purpose:   create selected_events.txt file out of pdata_events.txt 
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
from obspy.imaging.beachball import Beach
import os

# ------------------- INPUT -----------------------------
min_mag = 6.0
max_mag = 9.999999999
plot_ev = True
# -------------------------------------------------------

fio_pdata_events = open(os.path.join('.', 'results', 'pdata_events.txt'), 'r')
pdata_events = fio_pdata_events.readlines()[1:]
print '%s total number of events' %(len(pdata_events))
fio_selected_events = open(os.path.join('.', 'results', 'selected_events.txt'), 'w')
enum = 0
for i in range(len(pdata_events)):
    if min_mag<=float(pdata_events[i].split(',')[4])<=max_mag:
        fio_selected_events.writelines(pdata_events[i])
        enum += 1

print '%s found events' %(enum)
fio_selected_events.close()

if plot_ev:
    print 'Plotting the events!'
    m = Basemap(projection='cyl', lon_0=0.0, lat_0=0.0, resolution='c')
    m.drawcoastlines()
    m.fillcontinents()
    m.drawparallels(np.arange(-90., 120., 30.))
    m.drawmeridians(np.arange(0., 420., 60.))
    m.drawmapboundary()
    fio_selected_events = open(os.path.join('.', 'results', 'selected_events.txt'), 'r')
    selected_events = fio_selected_events.readlines()
    for i in range(len(selected_events)):
        try:
            evnt = selected_events[i].split(',')
            x, y = m(float(evnt[2]), float(evnt[1]))
            focmecs = [float(evnt[5]),float(evnt[6]),float(evnt[7]),\
                        float(evnt[8]),float(evnt[9]),float(evnt[10])]
            ax = plt.gca()
            b = Beach(focmecs, xy=(x, y), width=3, linewidth=1, alpha=0.85)
            b.set_zorder(10)
            ax.add_collection(b)
        except Exception, e:
            print evnt[0], e

    plt.title('%s Events found for:\n%s<=magnitude<=%s' %(enum, min_mag, max_mag))
    plt.show()
