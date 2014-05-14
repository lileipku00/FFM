#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  yspec2grf.py
#   Purpose:   yspec2grf main program 
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python and Obspy modules will be imported in this part.
import glob
import numpy as np
from obspy.core import read, Trace, UTCDateTime
import pprocess
import subprocess
import sys, os
import shutil

'''
sys.argv[1]: path to the newly calculated YSPEC in the gallery
sys.argv[2]: number of processors 
sys.argv[3]: source time
sys.argv[4]: real station names file
sys.argv[5]: phase name
sys.argv[6]: data folder (where the results should move for FFM) 
sys.argv[7]: travel time obtained from MATLAB code
sys.argv[8]: solver
'''

t = UTCDateTime(sys.argv[3])

########################################################################
###################### Functions are defined here ######################
########################################################################
    
###################### epi_dist ###################################


def epi_dist(tr, req_phase='Pdiff', tb=20, ta=100, model='iasp91', forward_code='yspec'):
    """
    returns O, A, B, E, GCARC for SAC headers
    this will be done directly by reading A from SAC header!
    It should be completely compatible with WKBJ method
    """
    if not forward_code == 'yspec':
        print 'O, B and E headers should be set correctly!'
    if not tr.stats.sac.a == -12345.0:
        t_phase = tr.stats.sac.a
        phase_exist = 'Y'
        tr_sliced = tr.slice(tr.stats.starttime + t_phase - tb, tr.stats.starttime + t_phase + ta)
        # XXX O should be zero? for the moment it could be zero
        # since YSPEC seismograms are all start at 0 but not for AXISEM
        O = 0
        A = t_phase
        B = tr_sliced.stats.starttime - tr.stats.starttime
        # XXX Shouldn't be related to tr_sliced (end time)?
        E = tr_sliced.stats.endtime - tr.stats.starttime
        GCARC = tr.stats.sac.gcarc
    else:
        phase_exist = 'N'
        tr_sliced = 0 
        O = 0
        A = 0
        B = 0
        E = 0
        GCARC = 0
    
    return (phase_exist, O, A, B, E, GCARC, tr_sliced)

###################### y2m ###################################


def y2m(file, path):
    """
    yspec outputs to SAC format
    """
    stationID = int(file.split('.')[-1])
    chans = ['BHZ', 'BHN', 'BHE']
    dat = np.loadtxt(file)
    npts = len(dat[:,0])
    for i, chan in enumerate(chans):
        stats = {'network': 'SG',
                 'station': 'RS%02d' % stationID,
                 'location': '',
                 'channel': chan,
                 'npts': npts,
                 'sampling_rate': (npts - 1.)/(dat[-1,0] - dat[0,0]),
                 'starttime': t,
                 'mseed': {'dataquality': 'D'}}
        traces = Trace(data=dat[:,1+i], header=stats)
        traces.write(os.path.join(path, 'SAC', 'dis.%s.%s.%s' % (traces.stats.station, traces.stats.location,
                                                                 traces.stats.channel)), format='SAC')

# ----------------- COLAT ---------------------
def COLAT(X):
    return 1.57079632-np.arctan([.993277*np.tan(X)])

# ----------------- AZIDL ---------------------
def AZIDL(SLAT, SLON, ELAT, ELON):
    SLAT = float(SLAT)*np.pi/180.
    SLON = float(SLON)*np.pi/180.
    ELAT = float(ELAT)*np.pi/180.
    ELON = float(ELON)*np.pi/180.
    SCOLAT=COLAT(SLAT)
    ECOLAT=COLAT(ELAT)
    SC1=np.sin(SCOLAT)
    SC2=np.cos(SCOLAT)
    SC3=np.sin(SLON)
    SC4=np.cos(SLON)
    EC1=np.sin(ECOLAT)
    EC2=np.cos(ECOLAT)
    EC3=np.sin(ELON)
    EC4=np.cos(ELON)
    AE=EC1*EC4
    BE=EC1*EC3
    AZI1=(AE-SC3)**2+(BE+SC4)**2+EC2*EC2-2.
    AZI2=(AE-SC2*SC4)**2+(BE-SC2*SC3)**2+(EC2+SC1)**2-2.

    if AZI2 == 0.:
        AZIS=np.pi-np.sign([AZI1])*np.pi/2.
    else:
        AZIS=np.arctan(complex(AZI2, AZI1))
        CODELB=SC1*(AE*SC4+BE*SC3)+SC2*EC2

    DEL=np.arccos(CODELB)
    AS=SC1*SC4
    BS=SC1*SC3
    AZI1=(AS-EC3)**2+(BS+EC4)**2+SC2*SC2-2.
    AZI2=(AS-EC2*EC4)**2+(BS-EC2*EC3)**2+(SC2+EC1)**2-2.

    if AZI2 == 0.:
        AZIE=np.pi-np.sign([AZI1])*np.pi/2.
    else:
        AZIE=np.arctan(complex(AZI2, AZI1))

    SCOLAT=np.pi/2.-SLAT
    ECOLAT=np.pi/2.-ELAT
    SC1=np.sin(SCOLAT)
    SC2=np.cos(SCOLAT)
    EC1=np.sin(ECOLAT)
    EC2=np.cos(ECOLAT)
    AE=EC1*EC4
    BE=EC1*EC3
    DELBOD=np.arccos(SC1*(AE*SC4+BE*SC3)+SC2*EC2)
    COSDEL=np.cos(DELBOD)
    X1=((SC2+EC2)**2)/(1.+COSDEL)
    X2=((SC2-EC2)**2)/(1.-COSDEL)
    X=X1+X2
    Y=X1-X2
    COTDEL=1./np.tan(DELBOD)
    SINDEL=np.sin(DELBOD)
    DEL2=DELBOD*DELBOD
    A=64.*DELBOD+16.*DEL2*COTDEL
    D=48.*SINDEL+8.*DEL2/COSDEL
    B=-(D+D)
    E=30.*np.sin(DELBOD+DELBOD)
    C=-30.*DELBOD-8.*DEL2*COTDEL-E/2.
    DELS=6378.2064*(DELBOD-.000847518825*(X*DELBOD-3.*Y*SINDEL)+.0897860195E-6*(X*(A+C*X+D*Y)+Y*(B+E*Y)))

    DEL_degree = DEL[0]*180./np.pi
    return DEL_degree

###################### td_modify ###################################


def td_modify(trace, req_phase, bg_model='iasp91'):
    """
    Modify gcarc and phase arrival based on ellipsoidal Earth model
    """
    ellip_gc = AZIDL(trace.stats.sac.stla, trace.stats.sac.stlo, trace.stats.sac.evla, trace.stats.sac.evlo)
    taup_process = subprocess.Popen(['taup_time', '-mod', bg_model, '-time', '-h', str(trace.stats.sac.evdp),
                                     '-ph', req_phase, '-deg', str(ellip_gc)], stdout=subprocess.PIPE)
    tt_raw = taup_process.communicate()[0]
    tt = tt_raw.split('\n')[0].split()[-1]
    if tt:
        tt = float(tt)
        if abs(tt - trace.stats.sac.a) >= 5:
            print 'Difference between ellipsoidal and spherical arrival time exceeds 5sec in %s.%s.%s.%s, proceed!' \
                  % (trace.stats.network, trace.stats.station, trace.stats.location, trace.stats.channel)
        else:
            trace.stats.sac.a = tt
            trace.stats.sac.gcarc = ellip_gc
    #else:
    #    print '%s does not exist in %s.%s.%s.%s' % (req_phase, trace.stats.network, trace.stats.station,
    #                                                trace.stats.location, trace.stats.channel)
    return trace

########################################################################
############################# Main Program #############################
########################################################################

req_phase = sys.argv[5]
forward_code = sys.argv[8]

print '\nConvert YSPEC outputs to SAC files!'
path1 = sys.argv[1]
if not os.path.isdir(os.path.join(path1, 'SAC')):
    os.mkdir(os.path.join(path1, 'SAC'))
    all_files = glob.glob(os.path.join(path1, 'yspec.out.*'))
    parallel_results = pprocess.Map(limit=int(sys.argv[2]), reuse=1)
    parallel_job = parallel_results.manage(pprocess.MakeReusable(y2m))
    for i in range(len(all_files)):
        parallel_job(all_files[i], path1)
    parallel_results.finish()
else:
    print '\nThe directory is already there:'
    print os.path.join(path1, 'SAC')


print '\nAdjusting the SAC names (fill up the header) according to the real data!'
path2 = sys.argv[4]
path3 = sys.argv[7]
if not os.path.isdir(os.path.join(path1, 'SAC_realName')):
    os.mkdir(os.path.join(path1, 'SAC_realName'))
    sta_name_open = open(path2)
    sta_name = sta_name_open.readlines()
    for i in range(len(sta_name)):
        sta_name[i] = sta_name[i].split(',')
        for j in range(len(sta_name[i])):
            sta_name[i][j] = sta_name[i][j].strip()
    for i in range(len(sta_name)):
        for chan in ['BHE', 'BHN', 'BHZ']:
            tr = read(os.path.join(path1, 'SAC', 'dis.RS' + '%02d' % (i+1) + '..' + chan))[0]
            tr.write(os.path.join(path1, 'SAC_realName', 'grf.%s.%s.%s.x00.%s' % (sta_name[i][0], sta_name[i][1],
                                                                                  sta_name[i][2], chan)), format='SAC')
            tr_new = read(os.path.join(path1, 'SAC_realName', 'grf.%s.%s.%s.x00.%s' % (sta_name[i][0], sta_name[i][1],
                                                                                       sta_name[i][2], chan)))[0]
            tr_new.stats.network = sta_name[i][0]
            tr_new.stats.station = sta_name[i][1]
            tr_new.stats.location = sta_name[i][2]
            tr_new.stats.channel = chan
            tr_new.stats.sac.stla = float(sta_name[i][5])
            tr_new.stats.sac.stlo = float(sta_name[i][6])
            tr_new.stats.sac.stel = float(sta_name[i][7])
            tr_new.stats.sac.stdp = float(sta_name[i][8])
            
            tr_new.stats.sac.evla = float(sta_name[i][9])
            tr_new.stats.sac.evlo = float(sta_name[i][10])
            tr_new.stats.sac.evdp = float(sta_name[i][11])
            tr_new.write(os.path.join(path1, 'SAC_realName', 'grf.%s.%s.%s.x00.%s' % (sta_name[i][0], sta_name[i][1],
                                                                                      sta_name[i][2], chan)),
                         format='SAC')
else:
    print '\nThe directory is already there:'
    print os.path.join(path1, 'SAC_realName')

print '\nStart filling in the arrival time and gcarc!' 
fio_ttime = open(path3, 'r')
fi_ttime = fio_ttime.readlines()
for _i in xrange(len(fi_ttime)):
    fi_ttime[_i] = fi_ttime[_i].split(',')[:-1]
for _i in xrange(len(fi_ttime)):
    try:
        tr = read(os.path.join(path1, 'SAC_realName', 'grf.%s.%s.%s.x00.%s' % (fi_ttime[_i][0], fi_ttime[_i][1],
                                                                               fi_ttime[_i][2], fi_ttime[_i][3])))[0]
        tr.stats.sac.gcarc = float(fi_ttime[_i][4])
        if not fi_ttime[_i][5] == 'NaN':
            tr.stats.sac.a = float(fi_ttime[_i][5])
        else:
            tr.stats.sac.a = -12345.0
        tr = td_modify(tr, req_phase)
        tr.write(os.path.join(path1, 'SAC_realName', 'grf.%s.%s.%s.x00.%s' % (fi_ttime[_i][0], fi_ttime[_i][1],
                                                                              fi_ttime[_i][2], fi_ttime[_i][3])),
                 format='SAC')
    except Exception, e:
        print 'ERROR: %s' % e

print '\nCutting Time-Window around %s' % req_phase
print '\nWARNING: tb=20, ta=100 (hard coded!)'
if not os.path.isdir(os.path.join(path1, 'grf_cut')):
    os.mkdir(os.path.join(path1, 'grf_cut'))

all_files = glob.glob(os.path.join(path1, 'SAC_realName', '*.*.*'))
for i in range(len(all_files)):
    tr = read(os.path.join(all_files[i]))[0]
    (phase_flag, O, A, B, E, GCARC, tr_sliced) = epi_dist(tr, req_phase=req_phase, tb=20, ta=100,
                                                          forward_code=forward_code)
    if phase_flag == 'Y':
        tr_sliced.stats.sac.o = O
        #tr_sliced.stats.sac.a = A
        tr_sliced.stats.sac.b = B
        tr_sliced.stats.sac.e = E
        #tr_sliced.stats.sac.gcarc = GCARC
        tr_sliced.write(os.path.join(path1, 'grf_cut', 'grf.%s.%s.%s.x00.%s' % (tr_sliced.stats.network,
                                                                                tr_sliced.stats.station,
                                                                                tr_sliced.stats.location,
                                                                                tr_sliced.stats.channel)),
                        format='SAC')

print '\nMove the data to data folder!'
data_dest = sys.argv[6]
if req_phase in ['P', 'Pdiff']:
    phase_ls = glob.glob(os.path.join(path1, 'grf_cut', '*.BHZ'))
    for fi in phase_ls:
        shutil.copy(fi, data_dest)
        os.remove(fi)
elif req_phase in ['SH']:
    print '\nit just moves the BHE and BHN!'
    phase_ls = glob.glob(os.path.join(path1, 'grf_cut', '*.BHE'))
    phase_ls.append(glob.glob(os.path.join(path1, 'grf_cut', '*.BHN')))
    for fi in phase_ls:
        shutil.copy(fi, data_dest)
        os.remove(fi)








# TRASH:

#self.dist = locations2degrees(lat1 = self.stats.sac.evla, 
#                long1 = self.stats.sac.evlo, 
#                lat2 = self.stats.sac.stla, 
#                long2 = self.stats.sac.stlo)
#epi_km = gps2DistAzimuth(tr.stats.sac.evla, tr.stats.sac.evlo,
#            tr.stats.sac.stla, tr.stats.sac.stlo)[0]
#dist = kilometer2degrees(epi_km/1000., radius = 6370.997)
#taup_process = subprocess.Popen(['taup_time', '-mod', model, '-time', '-h', 
#            str(tr.stats.sac.evdp), '-ph', req_phase, '-deg', str(dist)], stdout=subprocess.PIPE)
##exit_taup_code = os.waitpid(taup_process.pid, 0)
#tt_raw = taup_process.communicate()[0]
#tt = tt_raw.split('\n')[0]
#print tt
#if tt:
#    phase_exist = 'Y'
#    t_phase = float(tt)
#    tr_sliced = tr.slice(tr.stats.starttime + t_phase - tb, 
#            tr.stats.starttime + t_phase + ta)
#    O = 0
#    A = t_phase
#    B = tr_sliced.stats.starttime - tr.stats.starttime
#    E = tr_sliced.stats.endtime - tr.stats.starttime
#    GCARC = dist
#else:
#    phase_exist = 'N'
#    tr_sliced = 0 
#    O = 0; A = 0; B = 0; E = 0; GCARC = 0

#tt = getTravelTimes(dist, tr.stats.sac.evdp, model=model)
#t_phase = -12345.0
#for tt_item in tt:
#    if tt_item['phase_name'] == req_phase:
#        t_phase = tt_item['time']
#        phase_exist = 'Y'
#        break
#if t_phase != -12345.0:
#    #t_before = t_phase - tb
#    #t_after = t_phase + ta
#    tr_sliced = tr.slice(tr.stats.starttime + t_phase - tb, 
#            tr.stats.starttime + t_phase + ta)
#    O = 0
#    A = t_phase
#    B = tr_sliced.stats.starttime - tr.stats.starttime
#    E = tr_sliced.stats.endtime - tr.stats.starttime
#    GCARC = dist
#else:
#    phase_exist = 'N'
#    tr_sliced = 0 
#    O = 0; A = 0; B = 0; E = 0; GCARC = 0