#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  plot_dT_dA.py
#   Purpose:   plot dT and dA againts dominant period for FFM
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import glob
import os
import shutil
import sys

import util_ffproc as uf 

# ------------------- INPUT -----------------------------
xcorr_limit = 0.85
all_processed = True
update_all = False
plot_scatter = False
plot_mean = True
remote_dir = '/import/neptun-helles/hosseini/FFM' 
# -------------------------------------------------------

intro = 20*'-'
intro += '\nUsage:'
intro += '\npython plot_dT_dA.py bands[1-3] address-for-FFM-event.\n'
intro += 20*'-'

bands = sys.argv[1]
bands = range(int(bands[0]), int(bands[-1])+1)
band_period = {'1': 30.0,'2': 21.2,'3': 15.0,'4': 10.6,'5': 7.5,
                '6': 5.3,'7': 3.7,'8': 2.7}
if not all_processed:
    evadd = sys.argv[2]
    evname = evadd.split('/')[-1]
    if not evname: evname = evadd.split('/')[-2]
    all_staev = uf.reader(evadd, bands, band_period)
    passed_staev = uf.filters(all_staev, bands, xcorr_limit=xcorr_limit)
    if plot_scatter: uf.ffpscatter(passed_staev)
    if plot_mean: 
        per, dt_mean, a_mean, flag = uf.stamean(passed_staev)
        uf.ffplotmean(per, dt_mean)
else:
    proc_ev_ls = glob.glob(os.path.join(remote_dir, '*.*.*.*'))
    print '%s processed events found!' %(len(proc_ev_ls))
    if update_all:
        shutil.rmtree(os.path.join('.', 'statistics'))
        os.mkdir(os.path.join('.', 'statistics'))
    all_passed_staev = []
    for i in range(len(proc_ev_ls)):
        all_staev = uf.reader(proc_ev_ls[i], bands, band_period)
        if all_staev == []: continue
        passed_staev = uf.filters(all_staev, bands, xcorr_limit=xcorr_limit)
        if passed_staev[0] == []: continue
        all_passed_staev.append(passed_staev)
        if update_all: uf.writer(passed_staev, bands)
    if plot_scatter: uf.ffpscatter(all_passed_staev, all_events=True)
    if plot_mean:
        all_dt_mean = []
        all_a_mean = []
        for i in range(len(all_passed_staev)):
            per, dt_mean, a_mean, flag = uf.stamean(all_passed_staev[i])
            if flag:
                all_dt_mean.append([dt_mean, len(all_passed_staev[i][0])])
                all_a_mean.append([a_mean, len(all_passed_staev[i][0])])
        #uf.allffplotmean(per, all_dt_mean)
        uf.meanall_ffplot(per, all_dt_mean, all_a_mean)

