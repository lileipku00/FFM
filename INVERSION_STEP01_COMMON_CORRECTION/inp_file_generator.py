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
import datetime

import output_reader as outread

"""
TODO
- characterization of noise?
- check all the parameters with: /home/hosseini/Work/Scripts/FFT/UTILS/FFINVERSION/AMPLITUDES/Programs/write_raydata_input
"""

# ------------------- INPUT -----------------------------
events_dir = '/home/hosseini/Work/Scripts/gitHUB/MEASUREMENTS/P_measure_1_sec_LAMBDA_1-5_32_100'
phase = 'P'
req_band = 'band01'
all_events = '0274.2009.273.a'
#all_events = True

# Input file for raydata:
input_file_name_part = 'global'
input_counter = 1
twinned = 'None'

############## CRITERIA ###############
# WARNING: all the intervals are min <= ... < max
min_depth = -10
max_depth = 1000

min_xcorr = 0.85
max_xcorr = 1

min_epi = 32
max_epi = 85.01

check_clip = True
#######################################

bg_model = 'IASP91'
selected_events_add = './info/selected_events_indexed.txt'
# -------------------------------------------------------
t_start = datetime.datetime.now()

input_file_name = '%s_%s_%s' % (input_file_name_part, req_band, input_counter)
passed_event_adds = outread.event_filter(events_dir, selected_events_add=selected_events_add, all_events=all_events,
                                         min_dp=min_depth, max_dp=max_depth)

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

passed_ev_stas_array = np.array([])
for i in range(len(output_files)):
    if passed_ev_stas_array.size == 0:
        passed_ev_stas_array = output_files[i]
    else:
        passed_ev_stas_array = np.append(passed_ev_stas_array, output_files[i], 0)

filt_array = outread.array_station_filter(passed_ev_stas_array, min_xcorr=min_xcorr, max_xcorr=max_xcorr,
                                          min_epi=min_epi, max_epi=max_epi, check_clip=check_clip)
outread.check_selection(filt_array, min_xcorr, max_xcorr, min_epi, max_epi)
outread.raydata_input_generator(filt_array=filt_array, input_file_name=input_file_name, twinned=twinned, phase=phase,
                                min_xcorr=min_xcorr, min_depth=min_depth, max_depth=max_depth, min_epi=min_epi,
                                max_epi=max_epi, check_clip=check_clip)
outread.raydata_input(bg_model=bg_model, input_file_name=input_file_name, phase=phase)

print '#event-station pairs: %s' % len(filt_array)
print 'REGULAR END --- %s sec' % (datetime.datetime.now() - t_start)




#### ------------------------------- TRASH ----------------------------
####passed_ev_stas = []
###passed_ev_stas_array = np.array([])
###for i in range(len(output_files)):
###    #passed_ev_sta_single = outread.station_filter(output_files[i], min_xcorr=min_xcorr, max_xcorr=max_xcorr,
###    #                                              min_epi=min_epi, max_epi=max_epi, check_clip=check_clip)
###    #if not passed_ev_sta_single.size == 0:
###    #    passed_ev_stas.append(passed_ev_sta_single)
###    #    if passed_ev_stas_array.size == 0:
###    #        passed_ev_stas_array = passed_ev_sta_single
###    #    else:
###    #        passed_ev_stas_array = np.append(passed_ev_stas_array, passed_ev_sta_single, 0)
###    if passed_ev_stas_array.size == 0:
###        passed_ev_stas_array = output_files[i]
###    else:
###        passed_ev_stas_array = np.append(passed_ev_stas_array, output_files[i], 0)