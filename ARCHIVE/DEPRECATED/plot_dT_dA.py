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
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# ------------------- INPUT -----------------------------
xcorr_limit = 0.9
plot_on = False
write_on = True
plt_all_stas = False

# ATTENTION:
# da in this script is A observed and not really da!
# -------------------------------------------------------

intro = 20*'-'
intro += '\nUsage:'
intro += '\npython plot_dT_dA.py bands[1-3] address-for-FFM-event.\n'
intro += 20*'-'
#print intro

bands = sys.argv[1]
evadd = sys.argv[2]
bands = range(int(bands[0]), int(bands[-1])+1)

band_period = {'1': 30.0,'2': 21.2,'3': 15.0,'4': 10.6,'5': 7.5,
                '6': 5.3,'7': 3.7,'8': 2.7}

evname = evadd.split('/')[-1]
if not evname:
    evname = evadd.split('/')[-2]
passed_staev = []
for i in bands:
    try:
        str_i = 'band0' + str(i)
        passed_staev_tmp = []
        all_dt_event = np.array([])
        all_da_event = np.array([])
        fio_dt = open(os.path.join(evadd, 'outfiles', 'ffproc.ampstt.' + str_i), 'r')
        f_dt = fio_dt.readlines()
        for j in range(2, len(f_dt)):
            info_dt = f_dt[j].split()
            xcorr = float(info_dt[2])
            da = float(info_dt[3])
            dt = float(info_dt[5])
            lat = float(info_dt[6])
            lon = float(info_dt[7])
            passed_staev_tmp.append([lat, lon, xcorr, band_period[str(i)]])
            all_dt_event = np.append(all_dt_event, dt)
            all_da_event = np.append(all_da_event, da)
        all_dt_median = all_dt_event - np.median(all_dt_event)
        all_da_median = all_da_event #- np.median(all_da_event)
        for k in range(len(all_dt_median)):
            passed_staev_tmp[k].insert(2, all_dt_median[k])
        for k in range(len(all_da_median)):
            passed_staev_tmp[k].insert(3, all_da_median[k])
        passed_staev.append(passed_staev_tmp)
    except Exception, e:
        print e

dt_all = []
da_all = []
m_dt_all = np.array([])
co_dt_all = np.array([])
m_da_all = np.array([])
co_da_all = np.array([])
for i in range(len(passed_staev[0])):
    flag = 'T'
    # required stations are those that can pass the xcorr test in all bands
    # therefore, it first checks the xcorr in all bands for one station
    # then continue calculating the RMS fit to the data
    for j in range(len(bands)):
        if passed_staev[j][i][4] < xcorr_limit:
            flag = 'F'
    if flag == 'F': continue
    x_sta = np.array([])
    dt_sta = np.array([])
    da_sta = np.array([])
    for j in range(len(bands)):
        x_sta = np.append(x_sta, passed_staev[j][i][5])
        dt_sta = np.append(dt_sta, passed_staev[j][i][2])
        # divided by 1.e9 because of the difference in scale (real and syn)
        da_sta = np.append(da_sta, passed_staev[j][i][3]/1.e9)
    
    # calculate the least-square solution for dT and dA in each station
    A = np.vstack([x_sta, np.ones(len(x_sta))]).T
    m_dt, co_dt = np.linalg.lstsq(A, dt_sta)[0]
    m_dt_all = np.append(m_dt_all, m_dt)
    co_dt_all = np.append(co_dt_all, co_dt)
    m_da, co_da = np.linalg.lstsq(A, da_sta)[0]
    m_da_all = np.append(m_da_all, m_da)
    co_da_all = np.append(co_da_all, co_da)
    
    dt_all.append(dt_sta)
    da_all.append(da_sta)

if write_on:
    fio = open(os.path.join('.', 'mean_values.txt'), 'a+')
    msg = '%s,%s,%s,%s,%s,%s,%s,%s\n' %(evname, np.mean(m_dt_all), np.mean(co_dt_all),
                                np.mean(m_da_all), np.mean(co_da_all), 
                                band_period[str(bands[0])], band_period[str(bands[-1])],
                                len(m_dt_all))
    fio.writelines(msg)

if plot_on:
    plt.ion()
    plt.subplot(2, 1, 1)
    for j in range(len(dt_all)):
        x_one_sta = []
        y_one_sta = []
        for i in range(len(bands)):
            x_one_sta.append(x_sta[i])
            y_one_sta.append(dt_all[j][i])
        if plt_all_stas:
            plt.plot(x_one_sta, y_one_sta, 'o-')
    plt.ylabel('Time difference (dT)', fontsize = 'x-large', weight = 'bold')
    plt.xticks(fontsize = 'x-large', weight = 'bold')
    plt.yticks(fontsize = 'x-large', weight = 'bold')
    pltitle = evname
    pltitle += '\nxcorr >= %s' %(xcorr_limit)
    plt.title(pltitle, fontsize = 'x-large', weight = 'bold')

    # plot mean value (dt)
    x_mean = (band_period[str(bands[0])], band_period[str(bands[-1])])
    y_mean = (np.mean(m_dt_all) * band_period[str(bands[0])] + np.mean(co_dt_all),
               np.mean(m_dt_all) * band_period[str(bands[-1])] + np.mean(co_dt_all)) 
    plt.plot(x_mean, y_mean, '--', linewidth = 3)

    plt.subplot(2, 1, 2)
    for j in range(len(da_all)):
        x_one_sta = []
        y_one_sta = []
        for i in range(len(bands)):
            x_one_sta.append(x_sta[i])
            y_one_sta.append(da_all[j][i])
        if plt_all_stas:
            plt.plot(x_one_sta, y_one_sta, 'o-')
    plt.ylabel('Observed Amplitude', fontsize = 'x-large', weight = 'bold')
    plt.xlabel('Dominant Period', fontsize = 'x-large', weight = 'bold')
    plt.xticks(fontsize = 'x-large', weight = 'bold')
    plt.yticks(fontsize = 'x-large', weight = 'bold')

    # plot mean value (da)
    x_mean = (band_period[str(bands[0])], band_period[str(bands[-1])])
    y_mean = (np.mean(m_da_all) * band_period[str(bands[0])] + np.mean(co_da_all),
               np.mean(m_da_all) * band_period[str(bands[-1])] + np.mean(co_da_all)) 
    plt.plot(x_mean, y_mean, '--', linewidth = 3)

    plt.show()
