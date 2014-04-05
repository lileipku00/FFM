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
# It should be changed to -100 (large negative number) or so for nr_cc!!
xcorr_limit = -100
remote_dir = '/import/neptun-helles/hosseini/FFM/Pdiff_measure_2_sec_LAMBDA_1-5_90_180'

nr_cc = True
line_plot = True

# All stations is already disabled, so the following flag does not change anything
all_stations = False
just_high_cc = xcorr_limit
remove_GSN_median = True
# -------------------------------------------------------


# =======================================================
# =================== FUNCTIONS =========================
# =======================================================

# ------------------- round_to --------------------------


def round_to(n, precission):
    """
    rounding the numbers!
    """
    correction = 0.5 if n >= 0 else -0.5
    rounded = int(n/precission+correction)*precission
    rounded2 = round(rounded, 6)
    return rounded2

# ------------------- nr_dt -----------------------------


def nr_dt(t_shift_array, max_ts=30., width=0.1, num_bands=1, enum=0, leg='default', line_plot=False):
    """
    histogram plot for all measured traveltime anomalies
    EXAMPLES:
    1) For nr_cc:
    max_ts=2., width=0.1
    * max_ts=2., width=0.01
    2) For nr_dt:
    max_ts=30., width=0.1
    max_ts=30., width=0.5
    """
    bins = np.arange(-int(max_ts), int(max_ts)+width, width)
    
    for i in range(len(bins)): 
        bins[i] = round(bins[i], 6)
    for i in range(len(t_shift_array)):
        t_shift_array[i] = round_to(t_shift_array[i], width)
    
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
        if not line_plot:
            if i == 0:
                plt.bar(left=bins[i]-width*(0.25*1.5-enum*1.5*0.5/num_bands), width=1.5*width/2./num_bands,
                        height=digit_count[str(i)], color=dic_color[str(enum)], edgecolor=dic_color[str(enum)],
                        label=leg)
            else:
                plt.bar(left=bins[i]-width*(0.25*1.5-enum*1.5*0.5/num_bands),  width=1.5*width/2./num_bands,
                        height=digit_count[str(i)], color=dic_color[str(enum)], edgecolor=dic_color[str(enum)])
        else:
            x_line.append(bins[i])
            y_line.append(digit_count[str(i)])
    if line_plot:
        plt.plot(x_line[9:-10], y_line[9:-10], lw=3.0, label=leg, color=dic_color[str(enum)])

    #TRASH:
    #if not nr_cc:
    #if not nr_cc:
    ### |----------|----------|
    ### -1         0          1
    ###     ---->     <----
    ###else:
    ###    for i in range(len(digit_list)):
    ###        if t_shift_array[i] >= 0.:
    ###            digit_list[i] = digit_list[i]-1
    
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
        # [bands[i]] is passed like this because reader gets list as an input
        # reader(all_stations=False, just_high_cc=False, remove_GSN_median=False):
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
        for k in range(len(all_passed_staev[j][0])):
            if not all_passed_staev[j][0][k] == []:
                if not nr_cc:
                    t_shift_array.append(all_passed_staev[j][0][k][2])
                else:
                    # keep the name as t_shift_array to not change the whole script!
                    # However, if nr_cc is selected, it will be number of stations vs 
                    # cross correlation coefficient
                    t_shift_array.append(all_passed_staev[j][0][k][4])
    print 'Length of all passed data: %s' % len(t_shift_array)
    
    import py2mat_mod
    py2mat_mod.py2mat(t_shift_array, 't_shift_array_%s' % bands[i], 't_shift_array_%s' % bands[i])
    
    nr_dt(t_shift_array, num_bands=len(bands), enum=i, leg=str(band_period[str(bands[i])]) + 's', line_plot=line_plot)

if nr_cc:
    #Pdiff
    #plt.vlines(x=0.8, ymin=0.0, ymax=80000, lw=2, linestyle='--')
    #plt.xlim(-1.1, 1.1)
    #plt.ylim(ymax=80000)
    plt.vlines(x=0.8, ymin=0.0, ymax=215000, lw=2, linestyle='--')
    plt.xlim(-1.1, 1.1)
    plt.ylim(ymax=10000)
    plt.xlabel('xcorrelation factor', fontsize = 'xx-large', weight = 'bold')
    plt.ylabel('nr of data', fontsize = 'xx-large', weight = 'bold')
    plt.xticks(np.arange(-1.0, 1.1, 0.2), fontsize = 'xx-large', weight = 'bold')
    plt.yticks(fontsize = 'xx-large', weight = 'bold')
    #plt.title(fontsize = 'xx-large', weight = 'bold')
    plt.legend(loc=2, prop={'size':22})
    plt.xlim(0.01, 1.01)
    plt.show()
else:
    plt.xlim(-3.75, 3.75)
    plt.xlabel('Time lag / s', fontsize = 'xx-large', weight = 'bold')
    plt.ylabel('nr of data', fontsize = 'xx-large', weight = 'bold')
    plt.xticks(np.arange(-3.0, 4.0, 1.0), fontsize = 'xx-large', weight = 'bold')
    plt.yticks(fontsize = 'xx-large', weight = 'bold')
#plt.title(fontsize = 'xx-large', weight = 'bold')
    plt.legend(prop={'size':22})
    plt.show()
