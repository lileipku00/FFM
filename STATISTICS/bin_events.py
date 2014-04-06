#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  bin_events.py
#   Purpose:   histogram plot out of pdata_events.txt 
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
min_mag = -3.0
max_mag = 9.999999999
# -------------------------------------------------------

# ------------------- round_to --------------------------
def round_to(n, precission):
    correction = 0.5 if n >= 0 else -0.5
    return int(n/precission+correction)*precission

# ------------------- nr_dt -----------------------------
def nr_mag(array, width=1., enum=0):
    '''
    histogram plot for events
    '''
    bins = np.arange(-10.0, 700.0+width, width)
    for i in range(len(array)):
        array[i] = round_to(array[i], width)
    digit = np.digitize(array, bins)
    digit_list = digit.tolist()
    
    # |----------|----------|
    # -1         0          1
    #     ---->     <----
    for i in range(len(digit_list)):
        if array[i] >= 0.:
            digit_list[i] = digit_list[i]-1

    digit_count = {}
    for i in range(0, len(bins)):
        digit_count[str(i)] = digit_list.count(i)
    
    #x_line = []
    #y_line = []
    #for i in range(0, len(bins)):
    #    x_line.append(bins[i])
    #    y_line.append(digit_count[str(i)])
    #plt.plot(x_line, y_line, lw=3.0, label=leg, color='blue')
    
    for i in range(0, len(bins)):
        plt.bar(left = bins[i]-width*(0.25*1.5-enum*1.5*0.5), 
                width = 1.5*width/2., height = digit_count[str(i)], 
                color = 'blue', edgecolor = 'blue')

#---------------------------------------------------------
#--------------------MAIN PROGRAM-------------------------
#---------------------------------------------------------
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

print 'Plotting the events!'
fio_selected_events = open(os.path.join('.', 'results', 'selected_events.txt'), 'r')
selected_events = fio_selected_events.readlines()
ev_depth = []
for i in range(len(selected_events)):
    try:
        evnt = selected_events[i].split(',')
        ev_depth.append(float(evnt[3]))
    except Exception, e:
        print evnt[0], e

print max(ev_depth)
nr_mag(ev_depth, width=1., enum=0)
plt.xlabel('Depth(km)', size='x-large', weight='bold')
plt.ylabel('Number of events', size='x-large', weight='bold')
plt.xticks(size='x-large', weight='bold')
plt.yticks(size='x-large', weight='bold')
plt.show()
