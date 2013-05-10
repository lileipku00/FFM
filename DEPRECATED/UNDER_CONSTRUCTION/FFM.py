#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  FFM.py
#   Purpose:   Finite Frequency Measurements
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#
#   Copyright (C) 2013 Kasra Hosseini
#-------------------------------------------------------------------

# Added this line for python 2.5 compatibility
from __future__ import with_statement
import ConfigParser
import glob
import matplotlib.pyplot as plt
import numpy as np
from obspy import read, UTCDateTime
from obspy.signal import xcorr
import os
import sys
import time

from util_FFM import * 
from quake_handler import quake_info

#----------------------------- input handler -----------------------------
config = ConfigParser.RawConfigParser()
inputpath = 'input.cfg'

class input_handler:
    def __init__(self, inputpath):
        self.inpath = inputpath
        self.config = config.read(os.path.join(os.getcwd(), self.inpath))
        self.refpath = eval(config.get('address', 'ref'))
        self.reflabel = config.get('address', 'ref_label')
        self.wf1path = eval(config.get('address', 'wf1'))
        self.wf1label = config.get('address', 'wf1_label')
        self.stfpath = eval(config.get('address', 'stf'))
        self.identity = eval(config.get('station', 'identity'))
        self.req_phase = eval(config.get('phase', 'req_phase'))
        self.tb = float(config.get('phase', 'tb'))
        self.ta = float(config.get('phase', 'ta'))
        self.model = eval(config.get('phase', 'model'))
        self.stfgauss = eval(config.get('preprocess', 'stf_gauss'))
        self.stf_halfdur = float(config.get('preprocess', 'stf_halfdur'))
        self.stfreal = eval(config.get('preprocess', 'stf_real'))
        self.resample = eval(config.get('preprocess', 'resample'))
        self.rtr = eval(config.get('preprocess', 'rtr'))
        self.rmean = eval(config.get('preprocess', 'rmean'))
        self.taper = eval(config.get('preprocess', 'taper'))
        self.hfreq = float(config.get('preprocess', 'hfreq'))
        self.lfreq = float(config.get('preprocess', 'lfreq'))
        self.usesac = eval(config.get('preprocess', 'usesac'))

#--------------------------- Main Program -------------------------
print '\n----------------------------------------'
print 'Finite Frequency Measurements (FFM)'
print '----------------------------------------\n'
print '----------------------------------------'
war='WARNING: quake file should be checked before starting the program\n'
war+='this is due to the creation of quake file!\n'
war+='it will be addressed in the next versions in a way that the info\n'
war+='will be read from the quake file and not waveforms!'
print war
print '----------------------------------------\n'
time.sleep(2)

plt.ion()
plt1 = plt.figure()
ax1 = plt1.add_subplot(111)
plt2 = plt.figure()
ax2 = plt2.add_subplot(111)
plt3 = plt.figure()
ax3 = plt3.add_subplot(111)

# create the input class
inp = input_handler(inputpath)

st_1_proc=[]
st_2_proc=[]
epidist=[]
est_arrival=[]
ref_path = glob.glob(os.path.join(inp.refpath, inp.identity))
if inp.reflabel == 'REAL':
    events, address_events = \
        quake_info(address=os.path.join(inp.refpath, os.path.pardir),
        target = 'info')
    print '* Event date-time: %s' %events[0]['datetime']
    dt_sr = events[0]['datetime'] - UTCDateTime(0)
else:
    dt_sr = 0

if inp.stfreal:
    print '* Use the stf saved in %s' %inp.stfpath
    STF_tr = readSTF(inp.stfpath)
    STF_tr.plot()
elif inp.stfgauss:
    print '* Use Gaussian STF'
        
print '-------------------'
print '* Number of stations: %s' %len(ref_path)
for add in xrange(len(ref_path)):
    sys.stdout.write('\r')
    sys.stdout.write("[%-100s] %d%%" % ('='*int(100.*(add+1)/len(ref_path)),
                                            100.*(add+1)/len(ref_path)))
    sys.stdout.flush()
    add_ref = ref_path[add]
    id_ref_tmp = add_ref.split('/')[-1].split('.')
    id_ref = id_ref_tmp[1] + '.' + id_ref_tmp[2] + '.' + id_ref_tmp[3]
    tr_ref = read(add_ref)[0]
    if inp.wf1path:
        try:
            tr_2 = read(os.path.join(inp.wf1path, '*.' + id_ref))[0]
        except:
            print 'Warning: could not read %s for wf1' %id_ref
            continue
    else:
        tr_2 = False
    if inp.stfgauss:
        sigma =  inp.stf_halfdur / np.sqrt(2.) / 3.5
        if not inp.reflabel == 'REAL':
            tr_ref = convSTF_gauss(tr=tr_ref, sigma=sigma)
        if tr_2:
            tr_2 = convSTF_gauss(tr=tr_2, sigma=sigma)
    elif inp.stfreal:
        if tr_2:
            tr_2 = convSTF(tr=tr_2, STF_tr=STF_tr)
    tr_1_proc = preprocessing(tr_ref, usesac = inp.usesac, tr_address=inp.refpath, 
                                resample=inp.resample, rtr=inp.rtr, rmean=inp.rmean, 
                                taper=inp.taper, hfreq=inp.hfreq, lfreq=inp.lfreq)
    if tr_2: tr_2_proc = preprocessing(tr_2, usesac = inp.usesac, tr_address=inp.refpath,
                                resample=inp.resample, rtr=inp.rtr, rmean=inp.rmean, 
                                taper=inp.taper, hfreq=inp.hfreq, lfreq=inp.lfreq)

    if inp.reflabel == 'REAL':
        tr_cut = time_window(tr_1_proc)
        cut_info = tr_cut.epi_dist(req_phase=inp.req_phase, tb=inp.tb, 
                                    ta=inp.ta, model=inp.model)
        if cut_info[0] == 'Y':
            tr_1_proc.trim(starttime=events[0]['datetime']+cut_info[2], \
                            endtime=events[0]['datetime']+cut_info[3])
            try: 
                tr_1_proc.normalize()
                st_1_proc.append(tr_1_proc)
                est_arrival.append(cut_info[1])
                if tr_2: 
                    tr_2_proc.trim(starttime=UTCDateTime(0)+cut_info[2], \
                                    endtime=UTCDateTime(0)+cut_info[3])
                    tr_2_proc.normalize()
                    st_2_proc.append(tr_2_proc)
                epidist.append(cut_info[4])
            except: pass
    else:
        st_1_proc.append(tr_1_proc)
        st_2_proc.append(tr_2_proc)
        epidist.append(add)

#if inp.reflabel == 'REAL':
#    for tr in st_1_proc:
#        tr.data *= 1.e-9
# find the maximum absolute value in all the traces
#maxi=0
#for i in xrange(len(st_1_proc)):
#    maxi=max(maxi, np.abs(st_1_proc[i]).max(), np.abs(st_2_proc[i]).max())
maxi=0.1

# find the minimum time for plotting
mini_time = st_1_proc[0].stats.starttime
for i in xrange(len(st_1_proc)):
    ##change name in axisem does not work properly...so first it should be debugged and then:
    #st_1_id = '%s.%s.%s' %(st_1_proc[i].stats.station, st_1_proc[i].stats.location, st_1_proc[i].stats.channel)
    #st_2_id = '%s.%s.%s' %(st_2_proc[i].stats.station, st_2_proc[i].stats.location, st_2_proc[i].stats.channel)
    #if not st_1_id == st_2_id:
    #    print 'ERROR: %s and %s are not the same!' %(st_1_id, st_2_id)
    #    sys.exit(2)
    if tr_2: 
        mini_time=min(mini_time, st_1_proc[i].stats.starttime, 
                            st_2_proc[i].stats.starttime + dt_sr)
    else:
        mini_time=min(mini_time, st_1_proc[i].stats.starttime)
    
print '\n* Number of waveforms to be plotted: %s' %len(st_1_proc)
print '-------------------'
cc_info = []
for i in xrange(len(st_1_proc)):
    dt=st_1_proc[i].stats.delta
    npts=st_1_proc[i].stats.npts
    #t_1 = np.linspace(0., dt*(npts-1), npts) + (st_1_proc[i].stats.starttime-mini_time)
    t_1 = np.linspace(0., dt*(npts-1), npts)
    ax1.plot(t_1, st_1_proc[i].data/maxi+epidist[i], 'b', label=inp.reflabel)
    ax1.text(t_1[-1] + 0.2, epidist[i], '%s.%s.%s.%s'
                        %(st_1_proc[i].stats.network, st_1_proc[i].stats.station,
                        st_1_proc[i].stats.location, st_1_proc[i].stats.channel)) 
    if tr_2:
        dt=st_2_proc[i].stats.delta
        npts=st_2_proc[i].stats.npts
        #t_2 = np.linspace(0., dt*(npts-1), npts) + (st_2_proc[i].stats.starttime + dt_sr-mini_time)
        t_2 = np.linspace(0., dt*(npts-1), npts)
        ax1.plot(t_2, st_2_proc[i].data/maxi+epidist[i], 'gray', label=inp.wf1label)
        #Cross Correlation Coefficient
        t_cross, coeff_cross = xcorr(st_1_proc[i], st_2_proc[i],
                                    int(5.0*st_1_proc[i].stats.sampling_rate))
        cc_info.append(['%s.%s.%s.%s'
                        %(st_1_proc[i].stats.network, st_1_proc[i].stats.station,
                        st_1_proc[i].stats.location, st_1_proc[i].stats.channel),
                        t_cross, coeff_cross])
        ax2.scatter(coeff_cross, epidist[i], color = 'blue')
        ax2.axvline(x = 0.95, color = 'green', ls = 'dashed')
        ax3.scatter(t_cross/st_1_proc[i].stats.sampling_rate, epidist[i], color = 'red')
ax1.vlines(inp.tb, epidist[0]-3., epidist[-1], 'red')
l1, = ax1.plot(st_1_proc[0].data+epidist[0], 'b'); l1.remove()
if tr_2: 
    labels=[inp.reflabel, inp.wf1label]
    l2, = ax1.plot(st_2_proc[0].data+epidist[0], 'gray'); l2.remove()
    ax1.legend([l1, l2], labels)
else:
    labels=[inp.reflabel]
    ax1.legend([l1], labels)
plt.show()   
