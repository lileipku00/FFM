#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  create_events_file.py
#   Purpose:   create pdata_events.txt file out of pdata directory
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import glob
import os

#----------------- INPUT -------------------------
# remote pdata_processed directory
pdata_add = '/import/neptun-radler/AmplitudeProjects/pdata_processed/psdata_events'
# Work with specific events, 
# if target_events = all ---> it considers everything
target_events = 'all' 
#target_events = [\
#'0161.2009.196.a',
#'0161.2009.217.a',
#'0171.2009.230.a',
#'0173.2009.065.a',
#'0173.2009.211.a',
#'0174.2009.078.a',
#'0177.2009.144.a',
#'0178.2009.018.a',
#'0190.2009.174.a',
#'0192.2009.022.a',
#'0196.2009.003.a',
#'0284.2009.261.a',
#'0312.2009.264.a',
#'0325.2009.240.a',
#'0370.2009.182.c',
#'0381.2009.096.a',
#'0403.2009.157.a',
#'0412.2009.167.a',
#'0429.2009.085.a',
#'0429.2009.132.a',
#'0431.2009.155.a',
#'0641.2009.065.a',
#'0681.2009.188.a',
#'0684.2009.001.a',
#'0685.2009.260.a',
#'0691.2009.075.a',
#'0696.2009.074.a',
#'0702.2009.094.a',
#'0703.2009.222.a',
#'0718.2009.302.a',
#'0731.2009.148.a',
#'0756.2009.189.a']
#-------------------------------------------------

evs_ls = glob.glob(os.path.join(pdata_add, '*.*.*.*'))
print '%s events found in the archive!' %(len(evs_ls))
print '\nProblematic events:'
first_line = '#eventID,lat,lon,depth,mag,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp\n'
ev_info = []
for i in range(len(evs_ls)):
    try:
        ev = evs_ls[i]
        ev_name = ev.split('/')[-1]
        if target_events != 'all':
            if ev_name not in target_events:
                print '%s not in the target events...continue!' %(ev_name)
                continue
        fio_source = open(os.path.join(ev, 'outfiles', 'ampinv.source'), 'r')
        f_source = fio_source.readlines()
        ev_year, ev_julianday, ev_hr, ev_min, ev_sec, ev_msec = f_source[1].split()
        evlat, evlon, catalog_depth, inverted_depth = f_source[3].split()
        
        # if the source file does not contain the inverted source info,
        # the data will be read based on NEIC, HARVARD cataloges
        try:
            mrr, mtt, mpp, mrt, mrp, mtp = f_source[13].split()
            event_depth = inverted_depth
        except Exception, e:
            mrr, mtt, mpp, mrt, mrp, mtp = f_source[7].split()
            event_depth = catalog_depth
        fio_mag = open(os.path.join(ev, 'README'), 'r')
        f_mag = fio_mag.readlines()
        ev_mag = f_mag[1].split()[1]
        ev_info.append('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n' %(ev_name,evlat,
            evlon,event_depth,ev_mag,mrr,mtt,mpp,mrt,mrp,mtp))
    except Exception, e:
        print ev_name

ev_info.insert(0, first_line) 
if not os.path.isdir(os.path.join('.', 'results')):
    os.mkdir(os.path.join('.', 'results'))
fio_ls_event = open(os.path.join('.', 'results', 'pdata_events.txt'), 'w')
for i in ev_info:
    fio_ls_event.writelines(i)
fio_ls_event.close() 
