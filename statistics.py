#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  statistics.py
#   Purpose:   statistical analysis of FF measurements
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
import numpy as np
import os
import sys

import util_ffproc as uf 

# ------------------- INPUT -----------------------------
xcorr_limit = 0.85
remote_dir = '/import/neptun-radler/hosseini-downloads/KASRA/FFM' 
# -------------------------------------------------------

# ------------------- nr_dt -----------------------------
def nr_dt(t_shift_array, max_ts=14, width=1, num_bands=1, enum=0, leg='default'):
    '''
    histogram plot for all measured traveltime anomalies
    '''
    bins = np.arange(-int(max_ts), int(max_ts), width)
    digit = np.digitize(t_shift_array, bins)
    digit_list = digit.tolist()
    digit_count = {}
    for i in range(0, len(bins)):
        digit_count[str(i)] = digit_list.count(i)
    
    dic_color = {'0': 'blue', '1': 'deepskyblue', '2': 'green', 
                    '3': 'darkorange', '4': 'brown', 
                    '5': 'olive', '6': 'tan', '7': 'darkviolet'}
    for i in range(0, len(bins)):
        if i==0:
            plt.bar(left = bins[i]-width*(0.25-enum*0.5/num_bands), 
                    width = width/2./num_bands, height = digit_count[str(i)], 
                    color = dic_color[str(enum)], 
                    edgecolor = dic_color[str(enum)], 
                    label=leg)
        else:
            plt.bar(left = bins[i]-width*(0.25-enum*0.5/num_bands), 
                    width = width/2./num_bands, height = digit_count[str(i)], 
                    color = dic_color[str(enum)], 
                    edgecolor = dic_color[str(enum)]) 
    
# --------------------------------------------------------------
# ------------------- MAIN PROGRAM -----------------------------
# --------------------------------------------------------------
bands = sys.argv[1]
bands = range(int(bands[0]), int(bands[-1])+1)
band_period = {'1': 30.0,'2': 21.2,'3': 15.0,'4': 10.6,'5': 7.5,
                '6': 5.3,'7': 3.7,'8': 2.7}

proc_ev_ls = glob.glob(os.path.join(remote_dir, '*.*.*.*'))
print '%s processed events found!' %(len(proc_ev_ls))

for i in range(len(bands)):
    all_passed_staev = []
    for j in range(len(proc_ev_ls)):
        # [bands[i]] is defined like this because reader gets
        # list as an input
        all_staev = uf.reader(proc_ev_ls[j], [bands[i]], band_period)
        if all_staev == []: continue
        passed_staev = uf.filters(all_staev, [bands[i]], xcorr_limit=xcorr_limit)
        if passed_staev[0] == []: continue
        all_passed_staev.append(passed_staev)
    t_shift_array = []
    for j in range(len(all_passed_staev)):
        # [0] in all_passed_staev[j][0] shows the current band
        # in general we have a loop over bands and in each step there is just
        # one band that we are working with which is accessible by [0]
        for k in range(len(all_passed_staev[j][0])):
            if not all_passed_staev[j][0][k] == []:
                t_shift_array.append(all_passed_staev[j][0][k][2])
    nr_dt(t_shift_array, num_bands=len(bands), enum=i, leg=str(band_period[str(bands[i])]) + 's')

plt.xlabel('Time lag / s', fontsize = 'xx-large', weight = 'bold')
plt.ylabel('nr of data', fontsize = 'xx-large', weight = 'bold')
plt.xticks(fontsize = 'xx-large', weight = 'bold')
plt.yticks(fontsize = 'xx-large', weight = 'bold')
#plt.title(fontsize = 'xx-large', weight = 'bold')
plt.legend()
plt.show()
