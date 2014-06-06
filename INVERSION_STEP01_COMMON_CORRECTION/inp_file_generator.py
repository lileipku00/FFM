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
import multiprocessing
import numpy as np
import datetime

import output_reader as outread

"""
TODO
- characterization of noise?
- check all the parameters with: /home/hosseini/Work/Scripts/FFT/UTILS/FFINVERSION/AMPLITUDES/Programs/write_raydata_input
- Something is wrong with the station location
"""

# ------------------- INPUT -----------------------------
events_dir = '/home/hosseini/Work/Scripts/gitHUB/MEASUREMENTS/P_measure_1_sec_LAMBDA_1-5_32_100'
phase = 'P'
req_band = 'band01'
#all_events = ['0274.2009.273.a', '0173.2003.087.a', '0209.2005.016.a', '0306.2008.269.a', '0718.2004.072.a']
#all_events = ['0274.2009.273.a']
all_events = True

# Input file for raydata:
input_file_name_part = 'applied_ccorr_all'
twinned = 'None'

# Input file for raymatrix
vp_vs_Qs = [1, 0, 0]
kernel_quad_km = 20.0
vertex_file = 'vertices.USA10'
facet_file = 'facets.USA10'

# Parallel request:
parallel_exec = False
np_req = 4

############## CRITERIA ###############
# WARNING: all the intervals are min <= ... < max
min_depth = -10
max_depth = 1000

min_xcorr = 0.8
max_xcorr = 1.01

min_epi = 32
max_epi = 95.01

check_clip = True
#######################################

check_selections = True
run_raydata = True
run_raymatrix = False
corr_io_list = [1, 1, 1]     # Ellipticity, crustal correction, elevation
bg_model = 'IASP91'
selected_events_add = './info/selected_events_indexed.txt'
# -------------------------------------------------------
t_start = datetime.datetime.now()

print '\n======>> reading the event information and filter them %s <= depth < %s' % (min_depth, max_depth)
passed_event_adds = outread.event_filter(events_dir, selected_events_add=selected_events_add, all_events=all_events,
                                         min_dp=min_depth, max_dp=max_depth)

print '\n======>> reading ALL FFM outputs'

output_files = []
passed_ev_stas_array = np.array([])
for i in range(len(passed_event_adds)):
    output_file = outread.output_file_reader(evinfo=passed_event_adds[i], req_band=req_band)
    if output_file.size == 0:
        continue
    output_files.append(output_file)

    if passed_ev_stas_array.size == 0:
        passed_ev_stas_array = output_file
    else:
        passed_ev_stas_array = np.append(passed_ev_stas_array, output_file, 0)

print '\n========================='
print 'All events   : %s' % len(passed_event_adds)
print 'Used events  : %s' % len(output_files)
print 'Missed events: %s' % (len(passed_event_adds) - len(output_files))
print '========================='

print '\n======>> filtering FFM outputs:'
print '%s <= xcorr < %s' % (min_xcorr, max_xcorr)
print '%s <= epicentral < %s' % (min_epi, max_epi)
print 'check clip: %s' % check_clip

filt_array = outread.array_station_filter(passed_ev_stas_array, min_xcorr=min_xcorr, max_xcorr=max_xcorr,
                                          min_epi=min_epi, max_epi=max_epi, check_clip=check_clip)


if check_selections:
    print '\n======>> check the selected event-station pairs'
    outread.check_selection(filt_array, min_xcorr, max_xcorr, min_epi, max_epi)

print '\n========================'
print '#event-station pairs: %s' % len(filt_array)
print '========================'

print '\n======>> input files (raydata and raymatrix) generator'
if not parallel_exec:
    outread.compile_raydata_raymatrix()
    input_counter = 1
    input_file_name = '%s_%s_%s' % (input_file_name_part, req_band, input_counter)
    outread.raydata_input_generator(filt_array=filt_array, input_file_name=input_file_name, twinned=twinned,
                                    phase=phase, min_xcorr=min_xcorr, min_depth=min_depth, max_depth=max_depth,
                                    min_epi=min_epi, max_epi=max_epi, check_clip=check_clip)
    outread.raydata_input(bg_model=bg_model, input_file_name=input_file_name, phase=phase)
    outread.raymatrix_input(vp_vs_Qs=vp_vs_Qs, kernel_quad_km=kernel_quad_km, vertex_file=vertex_file,
                            facet_file=facet_file, input_file_name=input_file_name)
    print '\n======>> prepare output directory at: ./RESULTS/%s_dir' % input_file_name
    outread.prepare_dir(input_file_name=input_file_name)
    outread.run_raydata_raymatrix(input_file_name=input_file_name, raydata=run_raydata, raymatrix=run_raymatrix)
    filt_array_corr = outread.raydata_ccorr_reader(filt_array=filt_array, input_file_name=input_file_name,
                                                   corr_io_list=corr_io_list)
    outread.raydata_ccorr_writer(filt_array_corr=filt_array_corr, events_dir=events_dir)
    if run_raymatrix:
        outread.mat2asc_run(input_file_name=input_file_name)
else:
    outread.compile_raydata_raymatrix()
    if not np.mod(len(filt_array), np_req) == 0:
        sub_filt_array = np.split(filt_array[0:-np.mod(len(filt_array), np_req)], np_req)
        sub_filt_array.append(filt_array[-np.mod(len(filt_array), np_req):])
    else:
        sub_filt_array = np.split(filt_array, float(np_req))
    par_job = []
    for nj in range(len(sub_filt_array)):
        input_file_name = '%s_%s_%s' % (input_file_name_part, req_band, nj+1)
        par_job.append(multiprocessing.Process(target=outread.parallel_raydata_raymatrix,
                                               args=(sub_filt_array[nj], input_file_name, twinned, phase, min_xcorr,
                                                min_depth, max_depth, min_epi, max_epi, check_clip, bg_model,
                                                vp_vs_Qs, kernel_quad_km, vertex_file, facet_file)))
    for nj in range(len(par_job)):
        par_job[nj].start()
    outread.check_par_jobs(par_job, 0.2)

    for nj in range(len(sub_filt_array)):
        input_file_name = '%s_%s_%s' % (input_file_name_part, req_band, nj+1)
        outread.mat2asc_run(input_file_name=input_file_name)
    outread.vtk_generator(input_file_name_part=input_file_name_part, req_band=req_band, vertex_file=vertex_file,
                          facet_file=facet_file, parallel_exec=parallel_exec, len_dirs=len(sub_filt_array))


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