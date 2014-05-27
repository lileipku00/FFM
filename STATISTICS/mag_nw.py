#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  mag_nw.py
#   Purpose:   number of src-rcv pairs vs magnitude
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
xcorr_limit = 0.8
#remote_dir = '/home/hosseini/Work/Scripts/gitHUB/MEASUREMENTS/Pdiff_measure_1_sec_LAMBDA_1-5_90_180'
remote_dir = '/import/neptun-helles/hosseini/FFM_RESULTS/Pdiff_measure_1_sec_LAMBDA_1-5_90_180'

min_epi = 120
max_epi = 180
# -------------------------------------------------------
  
# ------------------- mag_finder -----------------------------


def mag_finder(ev_address, ev_fi):
    """
    find the magnitude of the event
    """
    ev_name = os.path.basename(ev_address)
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
band_period = {'1': 30.0, '2': 21.2, '3': 15.0, '4': 10.6, '5': 7.5, '6': 5.3, '7': 3.7, '8': 2.7}

proc_ev_ls = glob.glob(os.path.join(remote_dir, '*.*.*.*'))
print '%s processed events found!' % len(proc_ev_ls)

ev_fio = open('selected_events_ALL.txt')
ev_fi = ev_fio.readlines()
for i in range(len(ev_fi)):
    ev_fi[i] = [ev_fi[i].split(',')[0], ev_fi[i].split(',')[4]]

mag_dic = {}
for i in np.arange(5.0, 10.0, 0.5):
    mag_dic[str(i)] = [0, 0]

for i in range(len(bands)):
    all_passed_staev = []
    mag_all = np.array([])
    nw_all = np.array([])
    for j in range(len(proc_ev_ls)):
        # [bands[i]] is defined like this because reader gets
        # list as an input
        all_staev = uf.reader(proc_ev_ls[j], [bands[i]], band_period, all_stations=False, just_high_cc=False,
                              remove_GSN_median=False)
        if not all_staev:
            continue

        passed_staev = uf.filters(all_staev, [bands[i]], xcorr_limit=xcorr_limit)
        if not passed_staev[0]:
            continue
        magni = mag_finder(proc_ev_ls[j], ev_fi)
        mag_all = np.append(mag_all, magni)

        ### FIND EACH EVENT WITH SPECIFIED MAGNITUDE AND THE NUMBER OF WAVEFORNS ASSOCIATED ???????
        passed_staev_epi = []
        for k in range(len(passed_staev[0])):
            if min_epi <= passed_staev[0][k][6] < max_epi:
                passed_staev_epi.append(passed_staev[0][k])
        nw_all = np.append(nw_all, len(passed_staev_epi))
    for j in range(len(mag_all)):
        for k in mag_dic:
            if abs(uf.round_to(mag_all[j], 0.5) - float(k)) < 0.1:
                mag_dic[k][0] += nw_all[j]
                mag_dic[k][1] += 1
                break

    mag_sta_list = []
    for md in mag_dic:
        mag_sta_list.append([float(md), mag_dic[md][0], mag_dic[md][1]])
    mag_sta_list.sort()


    import py2mat_mod
    py2mat_mod.py2mat(mag_sta_list, 'mag_sta_%s' % bands[i], 'mag_sta_%s' % bands[i])

    plt.ion()
    plt.figure()
    plt.subplot(2, 1, 2)
    for j in mag_dic:
        if mag_dic[j][1] > 0.1:
            plt.bar(left=float(j)-0.05, width=0.1, height=mag_dic[j][0]/mag_dic[j][1])
            #plt.bar(left=float(j)-0.05, width=0.1, height=mag_dic[j][0])
    plt.xlabel('Magnitude', fontsize='xx-large', weight='bold')
    plt.ylabel('#waveforms/#events', fontsize='xx-large', weight='bold')
    plt.xticks(fontsize='xx-large', weight='bold')
    plt.yticks(fontsize='xx-large', weight='bold')
    plt.legend()
    plt.show()

    #plt.figure()
    plt.subplot(2, 1, 1)
    for j in mag_dic:
        if mag_dic[j][1] > 0.1:
            plt.bar(left=float(j)-0.05, width=0.1, height=mag_dic[j][0])
    plt.xlabel('Magnitude', fontsize='xx-large', weight='bold')
    plt.ylabel('#waveforms', fontsize='xx-large', weight='bold')
    plt.xticks(fontsize='xx-large', weight='bold')
    plt.yticks(fontsize='xx-large', weight='bold')
    plt.legend()
    plt.show()

print "Number of all events used: %s" % len(mag_all)

nr_all_sta = 0
for i in range(len(nw_all)):
    nr_all_sta += nw_all[i]
print "Number of all stations used: %s" % nr_all_sta
raw_input('Please press enter to exit the program!')
