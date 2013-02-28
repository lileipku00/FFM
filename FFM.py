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

#XXX gps2azimuth and find the available phases to cut!
#XXX cut the traces based on the phase time
#XXX apply cross correlation (simple one)

# Added this line for python 2.5 compatibility
from __future__ import with_statement
import ConfigParser
import glob
import matplotlib.pyplot as plt
import numpy as np
from obspy import read, UTCDateTime
import os
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
war='WARNING: quake file should be checked before starting the program\n'
war+='this is due to the creation of quake file!\n'
war+='it will be addressed in the next versions in a way that the info\n'
war+='will be read from the quake file and not waveforms!'
print war
print '----------------------------------------\n'
time.sleep(5)

# create the input class
inp = input_handler(inputpath)

st_1_proc = []
st_2_proc = []
ref_path = glob.glob(os.path.join(inp.refpath, inp.identity))
if inp.reflabel == 'REAL':
    events, address_events = \
        quake_info(address=os.path.join(inp.refpath, os.path.pardir),
        target = 'info')
    dt_sr = events[0]['datetime'] - UTCDateTime(0)
else:
    dt_sr = 0

if inp.stfreal:
    STF_tr = readSTF(inp.stfpath)
    STF_tr.plot()

for add in xrange(len(ref_path)):
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
    if inp.stfgauss and tr_2:
        print 'use Gaussian STF'
        sigma =  inp.stf_halfdur / np.sqrt(2.) / 3.5
        tr_ref = convSTF_gauss(tr=tr_ref, sigma=sigma)
        tr_2 = convSTF_gauss(tr=tr_2, sigma=sigma)
    elif inp.stfreal and tr_2:
        #STF_tr = readSTF(inp.stfpath)
        #STF_tr.plot()
        tr_2 = convSTF(tr=tr_2, STF_tr=STF_tr)
    tr_1_proc = preprocessing(tr_ref, usesac = inp.usesac, tr_address=inp.refpath, 
                                resample=inp.resample, rtr=inp.rtr, rmean=inp.rmean, 
                                taper=inp.taper, hfreq=inp.hfreq, lfreq=inp.lfreq)
    if tr_2: tr_2_proc = preprocessing(tr_2, usesac = inp.usesac, tr_address=inp.refpath,
                                resample=inp.resample, rtr=inp.rtr, rmean=inp.rmean, 
                                taper=inp.taper, hfreq=inp.hfreq, lfreq=inp.lfreq)
    st_1_proc.append(tr_1_proc)
    if tr_2: st_2_proc.append(tr_2_proc)


if inp.reflabel == 'REAL':
    for tr in st_1_proc:
        tr.data *= 1.e-9

# find the maximum absolute value in all the traces
#maxi=0
#for i in xrange(len(st_1_proc)):
#    maxi=max(maxi, np.abs(st_1_proc[i]).max(), np.abs(st_2_proc[i]).max())
maxi=1
[tr.normalize() for tr in st_1_proc]
if tr_2: [tr.normalize() for tr in st_2_proc]

# find the minimum time for plotting reasons
mini_time = st_1_proc[0].stats.starttime
for i in xrange(len(st_1_proc)):
    if tr_2: 
        mini_time=min(mini_time, st_1_proc[i].stats.starttime, st_2_proc[i].stats.starttime + dt_sr)
    else:
        mini_time=min(mini_time, st_1_proc[i].stats.starttime)

    
for i in xrange(len(st_1_proc)):
    dt=st_1_proc[i].stats.delta
    npts=st_1_proc[i].stats.npts
    t_1 = np.linspace(0., dt*(npts-1), npts) + (st_1_proc[i].stats.starttime-mini_time)
    plt.plot(t_1, st_1_proc[i].data/maxi+i, 'b', label=inp.reflabel)
    if tr_2:
        dt=st_2_proc[i].stats.delta
        npts=st_2_proc[i].stats.npts
        t_2 = np.linspace(0., dt*(npts-1), npts) + (st_2_proc[i].stats.starttime + dt_sr-mini_time)
        plt.plot(t_2, st_2_proc[i].data/maxi+i, 'gray', label=inp.wf1label)

l1, = plt.plot(st_1_proc[0].data, 'b'); l1.remove()
if tr_2: 
    labels=[inp.reflabel, inp.wf1label]
    l2, = plt.plot(st_2_proc[0].data, 'gray'); l2.remove()
    plt.legend([l1, l2], labels)
else:
    labels=[inp.reflabel]
    plt.legend([l1], labels)
plt.show()    
