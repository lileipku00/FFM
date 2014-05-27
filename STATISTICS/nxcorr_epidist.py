#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  nxcorr_epidist.py
#   Purpose:   plot xcorr vs epicentral distance of FF measurements
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
import pickle
import sys

import util_ffproc as uf 

# Calculating the percentage:
# Method: First we should run the code for xcorr_limit = -1000
# This can assure us that we consider all the source receiver pairs
# The next two lines should be uncomment first
# Second: continue after two lines....

#y_l = open(os.path.join('statistics_all', 'yband' + str(enum)), 'w')
#pickle.dump(y_line, y_l)

# (cont) in this part, first we load the data for all source receiver pairs...
# note that there is not difference between band01 to the last band! since we consider all
# Now uncomment the next lines (up to plt.plot)
# to load and calculate the percentage
# IMPORTANT: the two lines up should be commented out!

# ------------------- INPUT -----------------------------
xcorr_limit = -1000
#remote_dir = '/home/hosseini/Work/Scripts/gitHUB/MEASUREMENTS/Pdiff_measure_1_sec_LAMBDA_1-5_90_180'
remote_dir = '/import/neptun-helles/hosseini/FFM_RESULTS/Pdiff_measure_1_sec_LAMBDA_1-5_90_180'

# All stations is already disabled, so the following flag does not change anything
all_stations = False
just_high_cc = xcorr_limit
remove_GSN_median = True
# -------------------------------------------------------

# =======================================================
# ==================== FUNCTIONS ========================
# =======================================================

# ------------------- nr_dt -----------------------------


def nr_dt(t_shift_array, max_ts=180., width=1., num_bands=1, enum=0, leg='default', line_plot=False):
    """
    histogram plot for all measured traveltime anomalies
    """
    bins = np.arange(-int(max_ts), int(max_ts)+width, width)

    for i in range(len(bins)): 
        bins[i] = round(bins[i], 6)
    for i in range(len(t_shift_array)):
        t_shift_array[i] = uf.round_to(t_shift_array[i], width)
    
    digit = np.digitize(t_shift_array, bins)
    digit_list = digit.tolist()
   
    for i in range(len(digit_list)):
        digit_list[i] -= 1

    digit_count = {}
    for i in range(0, len(bins)):
        digit_count[str(i)] = digit_list.count(i)

    dic_color = {'0': 'blue', '1': 'deepskyblue', '2': 'green', '3': 'darkorange', '4': 'brown', '5': 'olive',
                 '6': 'tan', '7': 'darkviolet'}

    x_line = []
    y_line = []
    for i in range(0, len(bins)):
        x_line.append(bins[i])
        y_line.append(digit_count[str(i)])
    
    # Calculating the percentage:
    # Method: First we should run the code for xcorr_limit = -1000
    # This can assure us that we consider all the source receiver pairs
    # The next two lines should be uncomment first
    # Second: continue after two lines....
    
    y_l = open(os.path.join('statistics_all', 'yband' + str(enum)), 'w')
    pickle.dump(y_line, y_l)
    
    # (cont) in this part, first we load the data for all source receiver pairs...
    # note that there is not difference between band01 to the last band! since we consider all
    # Now uncomment the next lines (up to plt.plot)
    # to load and calculate the percentage
    # IMPORTANT: the two lines up should be commented out!
    
    #y_l = open(os.path.join('statistics_all', 'yband' + str(enum)))
    #y_l_all = pickle.load(y_l)
    #for i in range(len(y_l_all)-1, -1, -1):
    #    if y_l_all[i] < 1:
    #        del y_l_all[i]
    #        del y_line[i]
    #        del x_line[i]

    #plt.plot(x_line, np.array(y_line, dtype=float)/np.array(y_l_all, dtype=float)*100., lw=3.0, label=leg,
    #         color=dic_color[str(enum)])
    plt.plot(x_line, y_line, lw=3.0, label=leg, color=dic_color[str(enum)])
    
# --------------------------------------------------------------
# ------------------- MAIN PROGRAM -----------------------------
# --------------------------------------------------------------
bands = sys.argv[1]
bands = range(int(bands[0]), int(bands[-1])+1)
band_period = {'1': 30.0, '2': 21.2, '3': 15.0, '4': 10.6, '5': 7.5, '6': 5.3, '7': 3.7, '8': 2.7}

proc_ev_ls = glob.glob(os.path.join(remote_dir, '*.*.*.*'))
print '%s processed events found!' % len(proc_ev_ls)

for i in range(len(bands)):
    all_passed_staev = []
    for j in range(len(proc_ev_ls)):
        # [bands[i]] is defined like this because reader gets list as an input
        all_staev = uf.reader(proc_ev_ls[j], [bands[i]], band_period, all_stations=all_stations,
                              just_high_cc=just_high_cc, remove_GSN_median=remove_GSN_median)
        if not all_staev:
            continue
        passed_staev = uf.filters(all_staev, [bands[i]], xcorr_limit=xcorr_limit)
        if not passed_staev[0]:
            continue

        all_passed_staev.append(passed_staev)

    t_shift_array = []
    for j in range(len(all_passed_staev)):
        # [0] in all_passed_staev[j][0] shows the current band
        # in general we have a loop over bands and in each step there is just
        # one band that we are working with which is accessible by [0]
        # [6]: epicentral distance
        for k in range(len(all_passed_staev[j][0])):
            if not all_passed_staev[j][0][k] == []:
                t_shift_array.append(all_passed_staev[j][0][k][6])
    print 'Length of all passed data: %s' % len(t_shift_array)

    import py2mat_mod
    py2mat_mod.py2mat(t_shift_array, 't_shift_array_%s' % bands[i], 'xcorr_epi_t_shift_array_%s' % bands[i])

    nr_dt(t_shift_array, num_bands=len(bands), enum=i, leg=str(band_period[str(bands[i])]) + 's')

#plt.xlim(97.0, 160.0)
plt.xlim(30.0, 90.0)
plt.xlabel('Epicentral Distance / degree', fontsize = 'xx-large', weight = 'bold')
#plt.ylabel('nr of data', fontsize = 'xx-large', weight = 'bold')
plt.ylabel('% of data (xcorr>=0.8)', fontsize = 'xx-large', weight = 'bold')
plt.xticks(fontsize = 'xx-large', weight = 'bold')
plt.yticks(fontsize = 'xx-large', weight = 'bold')
#plt.title(fontsize = 'xx-large', weight = 'bold')
#plt.legend(prop={'size':22})
#plt.legend(prop={'size':22}, loc=8)
plt.show()



### #----------------------reader---------------------------------
### def reader(evadd, bands, band_period, all_stations=True, just_high_cc=False):
###     '''
###     This function reads the ffproc.ampstt.band....
###     '''
###     target_network = ['II', 'IU', 'CU', 'GT', 'IC']
###     passed_staev = []
###     for i in bands:
###         try:
###             str_i = 'band0' + str(i)
###             passed_staev_tmp = []
###             all_dt_event = np.array([])
###             all_da_event = np.array([])
###             all_dt_high_cc = []
###             fio_dt = open(os.path.join(evadd, 'outfiles',
###                             'ffproc.ampstt.' + str_i), 'r')
###             f_dt = fio_dt.readlines()
###             for j in range(2, len(f_dt)):
###                 info_dt = f_dt[j].split()
###                 if not all_stations:
###                     if not info_dt[9].split('.')[0] in target_network:
###                         continue
###                 xcorr = float(info_dt[2])
###                 da = float(info_dt[3])
###                 dt = float(info_dt[5])
###                 lat = float(info_dt[6])
###                 lon = float(info_dt[7])
###                 epi = float(info_dt[8])
###                 sta_id = info_dt[9]
###                 passed_staev_tmp.append([lat, lon, xcorr, band_period[str(i)], epi, sta_id])
###                 all_dt_event = np.append(all_dt_event, dt)
###                 all_da_event = np.append(all_da_event, da/1.e9)
###                 if just_high_cc:
###                     if xcorr >= just_high_cc:
###                         all_dt_high_cc.append(dt)
###             if just_high_cc and len(all_dt_high_cc) > 0:
###                 np_median = np.median(all_dt_high_cc)
###             else:
###                 np_median = np.median(all_dt_event)
###             all_dt_median = all_dt_event - np_median
###             all_da_median = all_da_event - np.median(all_da_event)
###             for k in range(len(all_dt_median)):
###                 passed_staev_tmp[k].insert(2, all_dt_median[k])
###             for k in range(len(all_da_median)):
###                 passed_staev_tmp[k].insert(3, all_da_median[k])
###             passed_staev.append(passed_staev_tmp)
###         except Exception, e:
###             print e
###     return passed_staev
