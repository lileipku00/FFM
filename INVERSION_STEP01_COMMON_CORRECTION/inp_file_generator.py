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

import output_reader as outread

# ------------------- INPUT -----------------------------
#events_dir = '/import/neptun-helles/hosseini/TEST/TEST_HERE'
events_dir = '/import/neptun-helles/hosseini/FFM_RESULTS/Pdiff_measure_1_sec_LAMBDA_1-5_90_180'
req_band = 'band01'

# WARNING: all the intervals are min <= ... < max
min_depth = 0
max_depth = 10

min_xcorr = 0.85
max_xcorr = 2

min_epi = 95
max_epi = 181

check_clip = True
# -------------------------------------------------------

passed_event_adds, passed_events_depths = outread.event_filter(events_dir, min_dp=min_depth, max_dp=max_depth)

output_files = []
for i in range(len(passed_event_adds)):
    output_file = outread.output_file_reader(add_output=passed_event_adds[i], req_band=req_band)
    if output_file.size == 0:
        continue
    output_files.append(output_file)

print '\n========================='
print 'All events   : %s' % len(passed_event_adds)
print 'Used events  : %s' % len(output_files)
print 'Missed events: %s' % (len(passed_event_adds) - len(output_files))
print '========================='

passed_ev_stas = []
for i in range(len(output_files)):
    passed_ev_sta_single = outread.station_filter(output_files[i], min_xcorr=min_xcorr, max_xcorr=max_xcorr,
                                                  min_epi=min_epi, max_epi=max_epi, check_clip=check_clip)
    if not passed_ev_sta_single.size == 0:
        passed_ev_stas.append(passed_ev_sta_single)
