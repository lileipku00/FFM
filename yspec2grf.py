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
from obspy.core import read, Trace, Stream, UTCDateTime
from obspy.core.util import gps2DistAzimuth
from obspy.core.util import kilometer2degrees
from obspy.taup.taup import getTravelTimes
import pprocess
import sys, os
import shutil
import subprocess

t = UTCDateTime(sys.argv[3])

########################################################################
###################### Functions are defined here ######################
########################################################################
    
###################### epi_dist ###################################
def epi_dist(tr, req_phase='Pdiff', tb=10, ta=25, model='iasp91'):
    
    if not tr.stats.sac.a == -12345.0:
        t_phase = tr.stats.sac.a
        phase_exist = 'Y'
        tr_sliced = tr.slice(tr.stats.starttime + t_phase - tb, 
                tr.stats.starttime + t_phase + ta)
        O = 0
        A = t_phase
        B = tr_sliced.stats.starttime - tr.stats.starttime
        E = tr_sliced.stats.endtime - tr.stats.starttime
        GCARC = tr.stats.sac.gcarc
    else:
        phase_exist = 'N'
        tr_sliced = 0 
        O = 0; A = 0; B = 0; E = 0; GCARC = 0

    
    
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

    
    return (phase_exist, O, A, B, E, GCARC, tr_sliced)


###################### y2m ###################################
def y2m(file, path):
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
                 'mseed' : {'dataquality': 'D'}}
        traces = Trace(data=dat[:,1+i], header=stats)
        traces.write(os.path.join(path, 'SAC', 'dis.%s.%s.%s' 
                        %(traces.stats.station, traces.stats.location, 
                        traces.stats.channel)), format = 'SAC')


########################################################################
############################# Main Program #############################
########################################################################
print '\nConvert YSPEC output to SAC files!'
path1 = sys.argv[1]
if not os.path.isdir(os.path.join(path1, 'SAC')):
    os.mkdir(os.path.join(path1, 'SAC'))
    all_files = glob.glob(os.path.join(path1, 'yspec.out.*'))
    parallel_results = pprocess.Map(limit=int(sys.argv[2]), reuse=1)
    parallel_job = parallel_results.manage(pprocess.MakeReusable(y2m))
    for i in range(0, len(all_files)):
        parallel_job(all_files[i], path1)
    parallel_results.finish()
else:
    print 'The directory is already there:'
    print os.path.join(path1, 'SAC')


print '\nAdjusting the SAC names (fill up the header) according to the real data!'
path2 = sys.argv[4]
path3 = sys.argv[7]
if not os.path.isdir(os.path.join(path1, 'SAC_realName')):
    os.mkdir(os.path.join(path1, 'SAC_realName'))
    sta_name_open = open(path2)
    sta_name = sta_name_open.readlines()
    for i in range(0, len(sta_name)):
        sta_name[i] = sta_name[i].split(',')
        for j in range(0, len(sta_name[i])):
            sta_name[i][j] = sta_name[i][j].strip()
    for i in range(0, len(sta_name)):
        for chan in ['BHE', 'BHN', 'BHZ']:
            tr = read(os.path.join(path1, 'SAC', 'dis.RS' + '%02d' % (i+1) + '..' + chan))[0]
            tr.write(os.path.join(path1, 'SAC_realName', 'grf.' + sta_name[i][0] + '.' + \
                sta_name[i][1] + '.' + sta_name[i][2] + '.x00.' + chan), format='SAC')
            tr_new=read(os.path.join(path1, 'SAC_realName', 'grf.' + sta_name[i][0] + '.' + \
                sta_name[i][1] + '.' + sta_name[i][2] + '.x00.' + chan))[0]
            tr_new.stats.network=sta_name[i][0]
            tr_new.stats.station=sta_name[i][1]
            tr_new.stats.location=sta_name[i][2]
            tr_new.stats.channel=chan
            tr_new.stats.sac.stla=float(sta_name[i][5])
            tr_new.stats.sac.stlo=float(sta_name[i][6])
            tr_new.stats.sac.stel=float(sta_name[i][7])
            tr_new.stats.sac.stdp=float(sta_name[i][8])
            
            tr_new.stats.sac.evla=float(sta_name[i][9])
            tr_new.stats.sac.evlo=float(sta_name[i][10])
            tr_new.stats.sac.evdp=float(sta_name[i][11])
            tr_new.write(os.path.join(path1, 'SAC_realName', 'grf.' + sta_name[i][0] + '.' + \
                sta_name[i][1] + '.' + sta_name[i][2] + '.x00.' + chan), format='SAC')
else:
    print 'The directory is already there:'
    print os.path.join(path1, 'SAC_realName')

print 'Start filling in the arrival time and gcarc!' 
fio_ttime = open(path3, 'r')
fi_ttime = fio_ttime.readlines()
for _i in xrange(len(fi_ttime)):
    fi_ttime[_i] = fi_ttime[_i].split(',')[:-1]
for _i in xrange(len(fi_ttime)):
    try:
        tr=read(os.path.join(path1, 'SAC_realName', 'grf.' + fi_ttime[_i][0] + '.' +
                fi_ttime[_i][1] + '.' + fi_ttime[_i][2] + '.x00.' + fi_ttime[_i][3]))[0]
        tr.stats.sac.gcarc = float(fi_ttime[_i][4])
        if not fi_ttime[_i][5] == 'NaN':
            tr.stats.sac.a = float(fi_ttime[_i][5])
        else:
            tr.stats.sac.a = -12345.0
        tr.write(os.path.join(path1, 'SAC_realName', 'grf.' + fi_ttime[_i][0] + '.' +
                fi_ttime[_i][1] + '.' + fi_ttime[_i][2] + '.x00.' + fi_ttime[_i][3]), 
                format='SAC')
    except Exception, e:
        print e

req_phase = sys.argv[5]
print '\nCutting Time-Window aournd %s' %req_phase
if not os.path.isdir(os.path.join(path1, 'grf_cut')):
    os.mkdir(os.path.join(path1, 'grf_cut'))

all_files = glob.glob(os.path.join(path1, 'SAC_realName', '*.*.*'))
for i in range(0, len(all_files)):
    tr = read(os.path.join(all_files[i]))[0]
    (phase_flag, O, A, B, E, GCARC, tr_sliced) = epi_dist(tr, req_phase=req_phase, tb=10, ta=25.6)
    if phase_flag == 'Y':
        tr_sliced.stats.sac.o = O
        #tr_sliced.stats.sac.a = A
        tr_sliced.stats.sac.b = B
        tr_sliced.stats.sac.e = E
        #tr_sliced.stats.sac.gcarc = GCARC
        tr_sliced.write(os.path.join(path1, 'grf_cut', 'grf.' + tr_sliced.stats.network + '.' +
            tr_sliced.stats.station + '.' + tr_sliced.stats.location + '.x00.' + 
            tr_sliced.stats.channel), format='SAC')

print '\nMove the data to data folder!'
data_dest = sys.argv[6]
shutil.rmtree(data_dest)
os.mkdir(data_dest)
if req_phase in ['P', 'Pdiff']:
    phase_ls = glob.glob(os.path.join(path1, 'grf_cut', '*.BHZ'))
    for fi in phase_ls:
        shutil.move(fi, data_dest)
elif req_phase in ['SH']:
    print 'it just move the BHE and BHN!'
    phase_ls = glob.glob(os.path.join(path1, 'grf_cut', '*.BHE'))
    phase_ls.append(glob.glob(os.path.join(path1, 'grf_cut', '*.BHN')))
    for fi in phase_ls:
        shutil.move(fi, data_dest)
