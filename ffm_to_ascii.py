#!/usr/bin/env python
# -*- coding: utf-8 -*-

#### XXXX IT HAS A READER: BE CAREFUL ABOUT XCORR AND MEDIAN!!!!

#-------------------------------------------------------------------
#   Filename:  ffm_to_ascii.py
#   Purpose:   Converts the FFM to ascii format (Ludwig)
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import glob
import numpy as np
import os
import sys

'''
FORMAT: (Ludwig)
ev_lon ev_lat depth sta_lon sta_lat dT
'''

# ------------------- INPUT -----------------------------
processed_events_add = '/import/neptun-helles/hosseini/FFM/P_measure_2_sec_LAMBDA_1-5'
xcorr_limit = 0.8

#----------------------reader---------------------------------
def reader(evadd, bands, band_period):
    '''
    This function reads the ffproc.ampstt.band....
    '''
    fio_source = open(os.path.join(evadd, 'infiles', 'yspec', 'in.yspec.source'))
    fi_source = fio_source.readlines()
    fio_source.close()
    evla, evlo, evdp = fi_source[0].split('\n')[0].split()
    #mrr, mtt, mpp, mrt, mrp, mtp = fi_source[1].split('\n')[0].split()
    #evyear, evmon, evday, evhour, evmin, evsec = \
    #    fi_source[3].split('\n')[0].split()
    
    passed_staev = []
    for i in bands:
        try:
            str_i = 'band0' + str(i)
            passed_staev_tmp = []
            all_dt_event = np.array([])
            all_da_event = np.array([])
            fio_dt = open(os.path.join(evadd, 'outfiles', 
                            'ffproc.ampstt.' + str_i), 'r')
            f_dt = fio_dt.readlines()
            for j in range(2, len(f_dt)):
                info_dt = f_dt[j].split()
                xcorr = float(info_dt[2])
                da = float(info_dt[3])
                dt = float(info_dt[5])
                lat = float(info_dt[6])
                lon = float(info_dt[7])
                sta_id = info_dt[9]
                passed_staev_tmp.append([lat, lon, xcorr, band_period[str(i)], sta_id])
                all_dt_event = np.append(all_dt_event, dt)
                all_da_event = np.append(all_da_event, da/1.e9)
            
            #for k in range(len(all_dt_event)):
            #    passed_staev_tmp[k].insert(2, all_dt_event[k])
            #    passed_staev_tmp[k].append([evla, evlo, evdp])
            all_dt_median = all_dt_event - np.median(all_dt_event)
            #all_da_median = all_da_event - np.median(all_da_event)
            
            for k in range(len(all_dt_median)):
                passed_staev_tmp[k].insert(2, all_dt_median[k])
                passed_staev_tmp[k].append([evla, evlo, evdp])
            #for k in range(len(all_da_median)):
            #    passed_staev_tmp[k].insert(3, all_da_median[k])
            passed_staev.append(passed_staev_tmp)
        except Exception, e:
            print e
    return passed_staev

#----------------------filters---------------------------------
def filters(all_staev, bands, xcorr_limit = False):
    '''
    filters all station-event pairs based on priori
    '''
    ind_failed = []
    for i in range(len(all_staev[0])):
        flag = 'T'
        # required stations are those that can pass the xcorr test in all bands
        # therefore, it first checks the xcorr in all bands for one station
        # then continue calculating the RMS fit to the data
        for j in range(len(bands)):
            if all_staev[j][i][3] < xcorr_limit:
                    flag = 'F'
            #if not all_stations:
            #    if not all_staev[j][i][6].split('.')[0] in target_network:
            #            flag = 'F'
        if flag == 'F':
            ind_failed.append(i)
    ind_failed.reverse()
    for j in range(len(bands)):
        for i in ind_failed:
            all_staev[j].pop(i)
    return all_staev 

# -------------------------------------------------------------------
# -------------------------- MAIN PROGRAM ---------------------------
# -------------------------------------------------------------------
bands = sys.argv[1]
bands = range(int(bands[0]), int(bands[-1])+1)
band_period = {'1': 30.0,'2': 21.2,'3': 15.0,'4': 10.6,'5': 7.5,
                '6': 5.3,'7': 3.7,'8': 2.7}

proc_ev_ls = glob.glob(os.path.join(processed_events_add, '*.*.*.*'))

all_passed_staev = []
for i in range(len(proc_ev_ls)):
    try:
        all_staev  = reader(proc_ev_ls[i], bands, band_period)
        if all_staev == []: continue
        passed_staev = filters(all_staev, bands, xcorr_limit=xcorr_limit) 
        if passed_staev[0] == []: continue
        all_passed_staev.append(passed_staev)
    except Exception, e:
        print e

ascii_fio = open(os.path.join('ascii-%s-%s.txt' %(bands[0], bands[-1])), 'w')

all_passed_staev.sort()
for i in range(len(all_passed_staev)):
    for j in range(len(all_passed_staev[i])):
        for k in range(len(all_passed_staev[i][j])):
            src_sta_pair = all_passed_staev[i][j][k]
            ascii_fio.writelines('%s %s %s %s %s %s\n' %(src_sta_pair[6][1], src_sta_pair[6][0], \
                                        src_sta_pair[6][2], src_sta_pair[1],  \
                                        src_sta_pair[0], src_sta_pair[2]))
ascii_fio.close()
