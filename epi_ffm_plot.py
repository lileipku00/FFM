#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-----------------------------------------------------------------------
#   Filename:  epi_ffm_plot.py
#   Purpose:   plotting tools for Finite Frequency Measurements 
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-----------------------------------------------------------------------

'''
===========
Instruction
===========

* Functionality: 
- this script plots the Matched Filter and real data for comparison purposes.
- the waveforms are arranged by epicentral distances.
* Input:
- list of the required stations to be plotted (in the same format as FFM output)
- event address
'''

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python and Obspy modules will be imported in this part.
import glob
import matplotlib.pyplot as plt
import numpy as np
from obspy import read, Trace, UTCDateTime
from obspy.taup import taup
from obspy.signal import xcorr
import os
import scipy.io
import sys

# ------------------- INPUT -----------------------------
# Inputs:
synthetic_add = '/import/neptun-helles/hosseini/FFM/YSPEC_SYN_GALLERY'
real_add = '/import/neptun-radler/AmplitudeProjects/psdata/' 
stf_add = '/import/neptun-radler/AmplitudeProjects/pdata_processed/psdata_events/'
# It is not important what we choose for ffproc_add, 
# since it just shows the stations available
ffproc_add = 'ffproc.ampstt.band01'
# normalize the data
normalize = True
# minimum and maximum epicentral distances for plotting
min_epi = -10.0
max_epi = 200.0
# phase
phase = 'P'
# time to cut the waveforms (tb: time before, ta: time after)
tb = 10
ta = 100
# apply a bandpass filter to the data
lfreq = 0.01
hfreq = 0.2
# minimum xcorrelation factor
min_xcorr = 0.85
# enum_max (how many waveforms should be considered)
enum_max = 100
# -------------------------------------------------------
# Try to read the event address
try: 
    ev_add = sys.argv[1]
except: sys.exit('python ffm_plot <address-to-data(FFM)>')

#--------------------------- preprocess --------------------------------
def preprocess(tr, lfreq, hfreq):
    '''
    all the pre-processing will be done here! 
    after this step the STF*tr will be calculated
    '''
    tr.filter('bandpass', freqmin=lfreq, freqmax=hfreq, corners=2, zerophase=True)
    return tr

#--------------------------- readSTF --------------------------------
def readSTF(add_stf, ev_time):
    """
    Read the Source Time Function
    Sampling rate is 10Hz
    """
    #reading and creating the Source Time Function
    stf_matlab = scipy.io.loadmat(add_stf)
    all_stf = stf_matlab['STF']
    all_time = stf_matlab['TAX']
    num_stf = len(all_stf[0])
    
    STF = []
    for j in range(num_stf):
        tmp_STF = np.array([all_stf[0][j]])
        for i in range(1, len(all_stf)):
            tmp_STF = np.append(tmp_STF, all_stf[i][j])
        stats = {'network': 'STF', 
                 'station': 'STF' + str(j), 
                 'location': '',
                 'channel': '00', 
                 'npts': len(all_time), 
                 'delta': (all_time[-1][j] - all_time[0][j])/(len(all_time)-1),
                 'starttime': (ev_time + all_time[0][j]),
                 'mseed' : {'dataquality': 'D'}}
        STF_tr = Trace(data=tmp_STF, header=stats)
        STF.append(STF_tr)
    return STF

#--------------------------- convSTF --------------------------------
def convSTF(tr, STF_tr, grp):
    """
    STF*Trace
    """
    # grp-1 is used here since group starts at 1!
    conv_tr = tr.copy()
    STF_tr[grp-1].resample(tr.stats.sampling_rate)
    conv_tr.data = np.convolve(STF_tr[grp-1].data, tr.data)
    return conv_tr

########################################################################
############################# Main Program #############################
########################################################################
# read the name of the stations from ffproc_add
sta_open = open(os.path.join(ev_add, 'outfiles', ffproc_add))
sta_read = sta_open.readlines()[2:]
for i in range(len(sta_read)):
    sta_read[i] = sta_read[i].split()

# determine the ev_name
ev_name = ev_add.split('/')[-1]
if not ev_name: ev_name = ev_add.split('/')[-2]

# temporarily read a green's fimctopm tp define the ev_time
# ev_time is required for constructing the Source Time Function (starttime)
grf_tmp = read(os.path.join(synthetic_add, ev_name, 'SAC_realName', 'grf.'+sta_read[0][-1]))[0]
ev_time = grf_tmp.stats.starttime
#print ev_time
STF_tr = readSTF(os.path.join(stf_add, ev_name, 'STF.mat'), ev_time)

# enum: enumerator, it controls how many waveforms should be plotted (enum_max)
enum = 0
sta_information = []
passed_nets = []
tr_all = []
for i in range(len(sta_read)):
    sta_name = sta_read[i][-1]
    grp = int(sta_read[i][1])
    
    # following parameters are not used yet in the program
    sta_lat = float(sta_read[i][-4])
    sta_lon = float(sta_read[i][-3])
    sta_epi = float(sta_read[i][-2])
    calc_xcorr = float(sta_read[i][2])
    calc_dt = float(sta_read[i][5])
    
    
    # read the green's function
    grf_tr = read(os.path.join(synthetic_add, ev_name, 'SAC_realName', 'grf.'+sta_name))[0]
    #event_1:
    #if not grf_tr.stats.station in ["109C", "AAK", "AAK", "ABKAR", "ACCN", "ADO", "AGMN", "AGRB", "AKGG", "AKUT", "ANMO", "ANWB", "BERG", "BESE", "BLO", "BLOW", "BMN", "BOAB", "BZN", "CCRK", "CMB", "DBO", "DLMT", "EUNU", "FACU", "FFD", "HDC", "HIA", "JSC", "LPAZ", "MGAN", "OBIP", "SC58", "SRU", "STD"]: continue
    #event_2:
    #if not grf_tr.stats.station in ["109C", "ACCN", "ACSO", "ADO", "ANMO", "ASCN", "ATE", "AUL", "BERG", "BESE", "BLO", "BLOW", "BMN", "BZN", "CCRK", "CMB", "CRAG", "DBIC", "DLMT", "DRLN", "DYA", "ESPZ", "FACU", "FFD", "GGNV", "LPW", "MTE", "PACT", "PESTR", "PPCWF", "PPTF", "RTC", "SC58", "SCO", "SRU", "STD"]: continue
    #if not grf_tr.stats.station in ["AAK", "AAK", "AKTO", "ANTO", "APE", "AQU", "HGN", "IPM", "KHC", "LSA", "MUN", "PSI", "SENIN", "SSB", "VNDA"]: continue
    #if grf_tr.stats.network in passed_nets: continue
    if min_epi <= grf_tr.stats.sac.gcarc <= max_epi:
        try:
            tt_list = taup.getTravelTimes(grf_tr.stats.sac.gcarc, grf_tr.stats.sac.evdp, model='iasp91')
            flag = 'searching'
            for j in range(len(tt_list)):
                if tt_list[j]['phase_name'] == phase:
                    phase_time = tt_list[j]['time']
                    #print '---------'
                    #print phase + ' is found and the arrival time:'
                    #print tt_list[j]['time']
                    #print '---------'
                    flag = 'found'
                    break
            if flag == 'searching': 
                print 'Could not find ' + phase + ' in: ' + str(grf_tr.stats.sac.gcarc)
                continue
            sta_name_split = sta_name.split('.')
            real_tr = read(os.path.join(real_add, ev_name, 'BH', 'dis.' + sta_name_split[1] + '.' + 
                                    sta_name_split[2] + '.' + sta_name_split[-1]))[0]
            mfi_tr = convSTF(grf_tr, STF_tr, grp)
            
            # preprocessing all three waveforms 
            mfi_tr = preprocess(tr=mfi_tr, lfreq=lfreq, hfreq=hfreq)
            real_tr = preprocess(tr=real_tr, lfreq=lfreq, hfreq=hfreq)
            real_tr.resample(mfi_tr.stats.sampling_rate) 
            
            real_tr = real_tr.slice(ev_time + phase_time - tb, ev_time + phase_time + ta)
            mfi_tr = mfi_tr.slice(ev_time + phase_time - tb, ev_time + phase_time + ta)
             
            np_xcorr, fac_xcorr = xcorr(mfi_tr.data, real_tr.data, int(15.*real_tr.stats.sampling_rate))
            
            sys.stdout.flush()
            #if fac_xcorr < min_xcorr: 
            #    print '.',
            #    continue
            #if np_xcorr == int(5.*real_tr.stats.sampling_rate): 
            #    print '.',
            #    continue
            
            print '\n******************'
            print sta_name
            print 'Sampling rate difference: ',
            print real_tr.stats.sampling_rate - mfi_tr.stats.sampling_rate
            print 'xcorr: ',
            print np_xcorr, fac_xcorr
            print 't phase: ',
            print phase_time
            print 'Enum: ',
            print enum 
            
            # following two lines control the maximum number of 'good' stations to be considered
            enum += 1
            if enum > enum_max: break
            t_diff = mfi_tr.stats.starttime - real_tr.stats.starttime
            
            #if normalize:
            #    t = np.linspace(0, (real_tr.stats.npts-1)/real_tr.stats.sampling_rate, real_tr.stats.npts)
            #    if i == 0: plt.plot(t, real_tr.data/abs(real_tr.max())+grf_tr.stats.sac.gcarc, 'black', label='REAL')
            #    else: plt.plot(t, real_tr.data/abs(real_tr.max())+grf_tr.stats.sac.gcarc, 'black')
            #    #plt.vlines(tb, grf_tr.stats.sac.gcarc-1, grf_tr.stats.sac.gcarc+1, linestyle='dashed')
            #    t = np.linspace(t_diff-np_xcorr/mfi_tr.stats.sampling_rate, t_diff - np_xcorr/mfi_tr.stats.sampling_rate + 
            #                        (mfi_tr.stats.npts-1)/mfi_tr.stats.sampling_rate, mfi_tr.stats.npts)
            #    #t = np.linspace(t_diff, t_diff + (mfi_tr.stats.npts-1)/mfi_tr.stats.sampling_rate, mfi_tr.stats.npts)
            #    #t = np.linspace(t_diff+calc_dt, t_diff + calc_dt + (mfi_tr.stats.npts-1)/mfi_tr.stats.sampling_rate, mfi_tr.stats.npts)
            #    if i == 0: plt.plot(t, mfi_tr.data/abs(mfi_tr.max())+grf_tr.stats.sac.gcarc, 'r', label='MFI')
            #    else: plt.plot(t, mfi_tr.data/abs(mfi_tr.max())+grf_tr.stats.sac.gcarc, 'r')
            tr_all.append([grf_tr.stats.sac.gcarc, real_tr, mfi_tr, t_diff, np_xcorr])
            sta_information.append([sta_name, sta_lat, sta_lon])
            passed_nets.append(grf_tr.stats.network)
        except Exception, e:
            print e

tr_all.sort(key=lambda x: int(x[0]))
for i in range(len(tr_all)):
    enum = i+1
    real_tr = tr_all[i][1]
    mfi_tr = tr_all[i][2]
    t_diff = tr_all[i][3]
    np_xcorr = tr_all[i][4]
    t = np.linspace(0, (real_tr.stats.npts-1)/real_tr.stats.sampling_rate, real_tr.stats.npts)
    if i == 0: plt.plot(t, real_tr.data/abs(real_tr.max()) + enum, 'black', label='REAL')
    else: plt.plot(t, real_tr.data/abs(real_tr.max()) + enum, 'black')
    #plt.vlines(tb, grf_tr.stats.sac.gcarc-1, grf_tr.stats.sac.gcarc+1, linestyle='dashed')
    t = np.linspace(t_diff-np_xcorr/mfi_tr.stats.sampling_rate, t_diff - np_xcorr/mfi_tr.stats.sampling_rate + 
                        (mfi_tr.stats.npts-1)/mfi_tr.stats.sampling_rate, mfi_tr.stats.npts)
    #t = np.linspace(t_diff, t_diff + (mfi_tr.stats.npts-1)/mfi_tr.stats.sampling_rate, mfi_tr.stats.npts)
    #t = np.linspace(t_diff+calc_dt, t_diff + calc_dt + (mfi_tr.stats.npts-1)/mfi_tr.stats.sampling_rate, mfi_tr.stats.npts)
    if i == 0: plt.plot(t, mfi_tr.data/abs(mfi_tr.max())+enum, 'r', label='MFI')
    else: plt.plot(t, mfi_tr.data/abs(mfi_tr.max())+enum, 'r')


plt.xlabel('Time (sec)') 
plt.ylabel('Epicentral Distance')
plt.ylim(ymax=enum+2)
plt.legend()

plt.show()
