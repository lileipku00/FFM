#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  mag_nw.py
#   Purpose:   number of waveforms vs magnitude
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
remote_dir = '/import/neptun-helles/hosseini/FFM' 
# -------------------------------------------------------
  
# ------------------- mag_finder -----------------------------
def mag_finder(ev_address, ev_fi):
    '''
    find the magnitude of the event
    '''
    ev_name = ev_address.split('/')[-1]
    for i in range(len(ev_fi)):
        if ev_name == ev_fi[i][0]:
            mag = float(ev_fi[i][1])
            break
    return mag

# --------------------------------------------------------------
# ------------------- MAIN PROGRAM -----------------------------
# --------------------------------------------------------------
bands = sys.argv[1]
bands = range(int(bands[0]), int(bands[-1])+1)
band_period = {'1': 30.0,'2': 21.2,'3': 15.0,'4': 10.6,'5': 7.5,
                '6': 5.3,'7': 3.7,'8': 2.7}

proc_ev_ls = glob.glob(os.path.join(remote_dir, '*.*.*.*'))
print '%s processed events found!' %(len(proc_ev_ls))

ev_fio = open('results/selected_events_all.txt')
ev_fi = ev_fio.readlines()
for i in range(len(ev_fi)):
    ev_fi[i] = [ev_fi[i].split(',')[0], ev_fi[i].split(',')[4]]

mag_dic = {}
for i in np.arange(5.0, 10.0, 0.1):
    mag_dic[str(i)] = [0, 0]

for i in range(len(bands)):
    all_passed_staev = []
    mag_all = np.array([])
    nw_all = np.array([])
    for j in range(len(proc_ev_ls)):
        # [bands[i]] is defined like this because reader gets
        # list as an input
        all_staev = uf.reader(proc_ev_ls[j], [bands[i]], band_period)
        if all_staev == []: continue
        passed_staev = uf.filters(all_staev, [bands[i]], xcorr_limit=xcorr_limit)
        if passed_staev[0] == []: continue
        magni = mag_finder(proc_ev_ls[j], ev_fi)
        mag_all = np.append(mag_all, magni)
        ### FIND EACH EVENT WITH SPECIFIED MAGNITUDE AND THE NUMBER OF WAVEFORNS ASSOCIATED ???????
        nw_all = np.append(nw_all, len(passed_staev[0]))
    for j in range(len(mag_all)):
        for k in mag_dic:
            if str(int(mag_all[j]*10.)/10.) == k:
                mag_dic[k][0] += nw_all[j]
                mag_dic[k][1] += 1
                break
    for j in mag_dic:
        if mag_dic[j][1] > 0.1:
            plt.bar(left=float(j)-0.05, width=0.1, height=mag_dic[j][0]/mag_dic[j][1])
            #plt.bar(left=float(j)-0.05, width=0.1, height=mag_dic[j][0])

#plt.xlim(-5.5, 5.5)
plt.xlabel('Magnitude', fontsize = 'xx-large', weight = 'bold')
plt.ylabel('#waveforms/#events', fontsize = 'xx-large', weight = 'bold')
plt.xticks(fontsize = 'xx-large', weight = 'bold')
plt.yticks(fontsize = 'xx-large', weight = 'bold')
#plt.title(fontsize = 'xx-large', weight = 'bold')
plt.legend()
plt.show()
