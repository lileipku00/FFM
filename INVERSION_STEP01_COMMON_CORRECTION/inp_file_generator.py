#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  inp_file_generator.py
#   Purpose:   Generate input file to be used by inversion step01 (common corrections)
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import numpy as np

import output_reader as outread

"""
TODO
- read event information from in.yspec.source
- event id: generate a list of all events and give an id to them, here use those ids
- kluster? I do not know how I can set this!
- characterization of noise?
- check all the parameters with: /home/hosseini/Work/Scripts/FFT/UTILS/FFINVERSION/AMPLITUDES/Programs/write_raydata_input
"""

# ------------------- INPUT -----------------------------
events_dir = '/home/hosseini/Work/Scripts/gitHUB/MEASUREMENTS/Pdiff_measure_1_sec_LAMBDA_1-5_90_180'
#events_dir = '/import/neptun-helles/hosseini/FFM_RESULTS/Pdiff_measure_1_sec_LAMBDA_1-5_90_180'
req_band = 'band01'

# WARNING: all the intervals are min <= ... < max
min_depth = -10
max_depth = 1000

min_xcorr = 0.8
max_xcorr = 2

min_epi = 0
max_epi = 181

check_clip = True
# -------------------------------------------------------

passed_event_adds = outread.event_filter(events_dir, min_dp=min_depth, max_dp=max_depth)

output_files = []
for i in range(len(passed_event_adds)):
    output_file = outread.output_file_reader(evinfo=passed_event_adds[i], req_band=req_band)
    if output_file.size == 0:
        continue
    output_files.append(output_file)

print '\n========================='
print 'All events   : %s' % len(passed_event_adds)
print 'Used events  : %s' % len(output_files)
print 'Missed events: %s' % (len(passed_event_adds) - len(output_files))
print '========================='

#passed_ev_stas = []
passed_ev_stas_array = np.array([])
for i in range(len(output_files)):
    #passed_ev_sta_single = outread.station_filter(output_files[i], min_xcorr=min_xcorr, max_xcorr=max_xcorr,
    #                                              min_epi=min_epi, max_epi=max_epi, check_clip=check_clip)
    #if not passed_ev_sta_single.size == 0:
    #    passed_ev_stas.append(passed_ev_sta_single)
    #    if passed_ev_stas_array.size == 0:
    #        passed_ev_stas_array = passed_ev_sta_single
    #    else:
    #        passed_ev_stas_array = np.append(passed_ev_stas_array, passed_ev_sta_single, 0)
    if passed_ev_stas_array.size == 0:
        passed_ev_stas_array = output_files[i]
    else:
        passed_ev_stas_array = np.append(passed_ev_stas_array, output_files[i], 0)
