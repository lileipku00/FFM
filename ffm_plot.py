#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-----------------------------------------------------------------------
#   Filename:  ffm_plot.py
#   Purpose:   plotting tools for Finite Frequency Measurements 
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python and Obspy modules will be imported in this part.
import glob
import matplotlib.pyplot as plt
import numpy as np
from obspy import read, Trace
import os
import scipy.io
import sys

# ------------------- INPUT -----------------------------
# Inputs:
gallery_add = '/import/neptun-helles/hosseini/FFM/YSPEC_SYN_GALLERY'
psdata_add = '/import/neptun-radler/AmplitudeProjects/psdata/' 
pdata_add = '/import/neptun-radler/AmplitudeProjects/pdata_processed/psdata_events/'
# scale the waveforms when it plots them
scale = 100
# normalize the data
normalize = True
# minimum and maximum epicentral distances for plotting
min_epi = 100.0
max_epi = 180.0
# time to cut the waveforms (tb: time before, ta: time after)
tb = 40
ta = 200
# apply a bandpass filter to the data
lfreq = 0.01
hfreq = 1.0
# -------------------------------------------------------

# Try to read the event address
try: ev_add = sys.argv[1]
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
def readSTF(add_stf):
    """
    Read the Source Time Function
    Sampling rate is 10Hz
    """
    #reading and creating the Source Time Function
    stf_matlab = scipy.io.loadmat(add_stf)
    stf = stf_matlab['STF']

    STF = np.array([stf[0][0]])
    for i in range(1, len(stf)):
        STF = np.append(STF, stf[i])
    
    stats = {'network': 'STF', 
             'station': 'STF', 
             'location': '',
             'channel': '01', 
             'npts': len(STF), 
             'sampling_rate': 10.0,
             'mseed' : {'dataquality': 'D'}}
    STF_tr = Trace(data=STF, header=stats)
    
    return STF_tr

#--------------------------- convSTF --------------------------------
def convSTF(tr, STF_tr):
    """
    STF*Trace
    """
    conv_tr = tr.copy()
    STF_tr.resample(tr.stats.sampling_rate)
    conv_tr.data = np.convolve(STF_tr.data, tr.data)
    return conv_tr

########################################################################
############################# Main Program #############################
########################################################################

# address for the green's functions and extracting the event name
grf_add = glob.glob(os.path.join(ev_add, 'data', '*.*.*'))
ev_name = ev_add.split('/')[-1]
if not ev_name:
    ev_name = ev_add.split('/')[-2]

# reading the STF (XXX: does not support several STF!)
STF_tr = readSTF(os.path.join(pdata_add, ev_name, 'STF.mat'))

if not normalize:
    maxi = 0
    print '\nStart to rescale the waveforms...',
    for add in grf_add:
        grf = read(add)[0]
        grf_name = add.split('/')[-1]
        tr_cmp = read(os.path.join(gallery_add, ev_name, 'SAC_realName', grf_name))[0]
        if min_epi <= grf.stats.sac.gcarc <= max_epi:
            maxi = max(abs(tr_cmp.max()), maxi)
    print 'DONE'
    print 'maximum: %s' %(maxi)

print 'length : %s' %(len(grf_add))
gca_arv = []

# grf: Green's function calculated and cut according to the arrival time
# tr_cmp: complete Green's function
# tr_real: real data
# tr_cmp_stf: complete Green's function convolved with STF
# grf_stf: Green's function convolved with STF
# XXX As a test here! the filter is applied once to tr_cmp and grf separately

for add in grf_add:
    print '.',
    sys.stdout.flush()
    grf = read(add)[0]
    print grf.stats.sac.gcarc
    if min_epi <= grf.stats.sac.gcarc <= max_epi:
        try:
            plt.clf()
            grf_name = add.split('/')[-1]
            tr_cmp = read(os.path.join(gallery_add, ev_name, 'SAC_realName', grf_name))[0]
            tr_cmp = tr_cmp.slice(grf.stats.starttime-tb, grf.stats.starttime+ta)
            tr_real = read(os.path.join(psdata_add, ev_name, 'BH', 'dis.%s.%s.%s'
                    %(tr_cmp.stats.station, tr_cmp.stats.location, tr_cmp.stats.channel)))[0]
            tr_real = tr_real.slice(tr_cmp.stats.starttime, tr_cmp.stats.endtime) 
            t_diff = grf.stats.starttime - tr_cmp.stats.starttime
            
            # preprocessing all three waveforms 
            grf = preprocess(tr=grf, lfreq=lfreq, hfreq=hfreq)
            tr_cmp = preprocess(tr=tr_cmp, lfreq=lfreq, hfreq=hfreq)
            tr_real = preprocess(tr=tr_real, lfreq=lfreq, hfreq=hfreq)
            
            # STF*GRF
            tr_cmp_stf = convSTF(tr_cmp, STF_tr)
            # After convolution with STF, it should be sliced again!?
            tr_cmp_stf = tr_cmp_stf.slice(grf.stats.starttime-tb, grf.stats.starttime+ta)
            grf_stf = convSTF(grf, STF_tr)
        except Exception, e:
            print e
        try:
            if normalize:
                #t = np.linspace(0, (tr_cmp.stats.npts-1)/tr_cmp.stats.sampling_rate, tr_cmp.stats.npts)
                #plt.plot(t, tr_cmp.data/abs(tr_cmp.max())+grf.stats.sac.gcarc, 'b', label='GRF', lw=2)
                t = np.linspace(0, (tr_real.stats.npts-1)/tr_real.stats.sampling_rate, tr_real.stats.npts)
                plt.plot(t, tr_real.data/abs(tr_real.max())+grf.stats.sac.gcarc, 'black', label='REAL', lw=2)
                #t = np.linspace(t_diff, (grf.stats.npts-1)/grf.stats.sampling_rate+t_diff, grf.stats.npts)
                #plt.plot(t, grf.data/abs(grf.max())+grf.stats.sac.gcarc, 'b', label='GRF(cut)')
                plt.vlines(tb + grf.stats.sac.a - grf.stats.sac.b, grf.stats.sac.gcarc-2, grf.stats.sac.gcarc+2, linestyle='dashed', lw=2)
                 
                t = np.linspace(0, (tr_cmp_stf.stats.npts-1)/tr_cmp_stf.stats.sampling_rate, tr_cmp_stf.stats.npts)
                plt.plot(t, tr_cmp_stf.data/abs(tr_cmp_stf.max())+grf.stats.sac.gcarc, 'r', label='STF*GRF', lw=2)
                
                #plt.xlim(0, 120)
                #plt.ylim(109.5, 111.5)
                #t = np.linspace(t_diff, (grf_stf.stats.npts-1)/grf_stf.stats.sampling_rate+t_diff, grf_stf.stats.npts)
                #plt.plot(t, grf_stf.data/abs(grf_stf.max())+grf_stf.stats.sac.gcarc, 'b', label='STF*GRF(cut)')

            else:
                t = np.linspace(0, (tr_cmp.stats.npts-1)/tr_cmp.stats.sampling_rate, tr_cmp.stats.npts)
                plt.plot(t, tr_cmp.data/maxi*scale+grf.stats.sac.gcarc, 'b', linestyle='dashed', label='GRF')
                t = np.linspace(0, (tr_real.stats.npts-1)/tr_real.stats.sampling_rate, tr_real.stats.npts)
                plt.plot(t, tr_real.data/maxi*scale/1.e9+grf.stats.sac.gcarc, 'black', label='REAL')
                t = np.linspace(t_diff, (grf.stats.npts-1)/grf.stats.sampling_rate+t_diff, grf.stats.npts)
                plt.plot(t, grf.data/maxi*scale+grf.stats.sac.gcarc, 'r', linestyle='dashed', label='GRF(cut)')
                plt.vlines(tb + grf.stats.sac.a - grf.stats.sac.b, grf.stats.sac.gcarc-20, grf.stats.sac.gcarc+20, linestyle='dashed')
                 
                t = np.linspace(0, (tr_cmp_stf.stats.npts-1)/tr_cmp_stf.stats.sampling_rate, tr_cmp_stf.stats.npts)
                plt.plot(t, tr_cmp_stf.data/maxi*scale+grf.stats.sac.gcarc, 'r', label='STF*GRF')
                t = np.linspace(t_diff, (grf_stf.stats.npts-1)/grf_stf.stats.sampling_rate+t_diff, grf_stf.stats.npts)
                plt.plot(t, grf_stf.data/maxi*scale+grf_stf.stats.sac.gcarc, 'b', label='STF*GRF(cut)')
        except Exception, e:
            print e
        
        #plt.xlim(0, 150)
        #plt.ylim(109.5, 111.5)
        plt.xlabel('Time (sec)', size='xx-large', weight='bold') 
        plt.ylabel('Epicentral Distance', size='xx-large', weight='bold')
        plt.xticks(size = 'large', weight = 'bold')
        plt.yticks(size = 'large', weight = 'bold')
        plt.title('%s\nEpicentral Distance: %s' %(grf_name, grf.stats.sac.gcarc))
        plt.legend()
        if not os.path.isdir(os.path.join('.', 'comparison_figs')):
            os.mkdir(os.path.join('.', 'comparison_figs'))
        plt.savefig(os.path.join('.', 'comparison_figs', grf_name + '.png'))
        plt.show()
        gca_arv.append([grf.stats.sac.gcarc, grf.stats.sac.a])

#plt.show()

#--------------------------- TRASH --------------------------------
#gca_arv.sort(key=lambda x: x[0])
#gca = []; arv_t = []
#for i in xrange(len(gca_arv)):
#    gca.append(gca_arv[i][0])
#    arv_t.append(gca_arv[i][1])

#plt.plot(arv_t, gca)

