#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  ffproc_reader.py
#   Purpose:   Read ffproc, process and generate the same file structure
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

import glob
import numpy as np
from obspy.core.util import gps2DistAzimuth
import os
import shutil
import sys

# ------------------- INPUT -----------------------------
remote_dir = './EVENTS'
output_dir = './ECORR'
phase = 'Pdiff'
single_event = False
# -------------------------------------------------------

# -------------------------- FUNCTIONS -----------------------------------


def reader_writer(evadd, bands, phase):
    """
    This function reads the ffproc.ampstt.band....
    """

    print '================================='
    print 'Start processing: %s' % evadd

    # meter to degrees converter
    m2degree = 360./(2.*np.pi*6371000.)

    failed = 0
    passed_staev = []
    for i in bands:
        try:
            str_i = 'band0' + str(i)
            passed_staev_band = []

            fio_source = open(os.path.join(evadd, 'outfiles', 'ampinv.source'), 'r')
            f_source = fio_source.readlines()
            fio_source.close()

            ev_year, ev_julianday, ev_hr, ev_min, ev_sec, ev_msec = f_source[1].split()
            evlat, evlon, catalog_depth, inverted_depth = f_source[3].split()
            try:
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[13].split()
            except Exception, e:
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[7].split()

            fio_dt = open(os.path.join(evadd, 'outfiles', 'ffproc.ampstt.' + str_i), 'r')
            f_dt = fio_dt.readlines()
            fio_dt.close()

            first_line = f_dt[0]
            second_line = f_dt[1]
            for j in range(2, len(f_dt)):
                info_dt = f_dt[j].split()
                indx = int(info_dt[0])
                grp = int(info_dt[1])
                xcorr = float(info_dt[2])
                da = float(info_dt[3])
                tobs = float(info_dt[4])
                dt = float(info_dt[5])
                stlat = float(info_dt[6])
                stlon = float(info_dt[7])
                epi = float(info_dt[8])
                sta_id = info_dt[9]
                evlat = float(evlat)
                evlon = float(evlon)
                gaz = gps2DistAzimuth(evlat, evlon, stlat, stlon)
                azimuth = gaz[1]
                epi_test = gaz[0]*m2degree
                assert(abs(epi - epi_test) < 1.0), 'WARNING: HUGE difference between read and calculated ' + \
                                                   'epicentral distance:\n%s - %s' % (epi, epi_test)
                passed_staev_band.append([indx, grp, xcorr, da, tobs, dt, stlat, stlon, epi, sta_id, evlat, evlon,
                                          inverted_depth, phase, azimuth])

            fio_ecorr_inp = open(os.path.join(os.path.curdir, '%s_band_%s.txt' % (os.path.basename(evadd), i)), 'w')

            for j in range(len(passed_staev_band)):
                passed_info = passed_staev_band[j]
                phase = passed_info[13]
                epi = float(passed_info[8])
                sdep = float(passed_info[12])
                slat = float(passed_info[10])
                azim = float(passed_info[14])
                tt = float(passed_info[5])

                if not j == len(passed_staev_band)-1:
                    fio_ecorr_inp.writelines('%s,%s,%s,%s,%s,%s\n' % (phase, epi, sdep, slat, azim, tt))
                else:
                    fio_ecorr_inp.writelines('%s,%s,%s,%s,%s,%s' % (phase, epi, sdep, slat, azim, tt))
            fio_ecorr_inp.close()

            fio_in_ecorr = open(os.path.join(os.path.curdir, 'in.ecorr'), 'w')
            fio_in_ecorr.writelines('%s_band_%s.txt\n' % (os.path.basename(evadd), i))
            fio_in_ecorr.writelines('%s_band_%s_output.txt' % (os.path.basename(evadd), i))
            fio_in_ecorr.close()

            os.system('./ecorr_blnk < in.ecorr')

            if not os.path.exists(os.path.join(output_dir, os.path.basename(evadd))):
                os.mkdir(os.path.join(output_dir, os.path.basename(evadd)))

            if not os.path.exists(os.path.join(output_dir, os.path.basename(evadd), 'outfiles')):
                os.mkdir(os.path.join(output_dir, os.path.basename(evadd), 'outfiles'))

            if not os.path.exists(os.path.join(output_dir, os.path.basename(evadd), 'outfiles', 'ampinv.source')):
                src = os.path.join(evadd, 'outfiles', 'ampinv.source')
                dst = os.path.join(output_dir, os.path.basename(evadd), 'outfiles', 'ampinv.source')
                shutil.copyfile(src, dst)

            if not os.path.exists(os.path.join(output_dir, os.path.basename(evadd), 'outfiles', 'ffproc.receivers')):
                src = os.path.join(evadd, 'outfiles', 'ffproc.receivers')
                dst = os.path.join(output_dir, os.path.basename(evadd), 'outfiles', 'ffproc.receivers')
                shutil.copyfile(src, dst)

            if not os.path.exists(os.path.join(output_dir, os.path.basename(evadd), 'outfiles', 'ffproc.source')):
                src = os.path.join(evadd, 'outfiles', 'ffproc.source')
                dst = os.path.join(output_dir, os.path.basename(evadd), 'outfiles', 'ffproc.source')
                shutil.copyfile(src, dst)

            fio_output = open(os.path.join(os.path.curdir, '%s_band_%s_output.txt' % (os.path.basename(evadd), i)), 'r')
            fi_output = fio_output.readlines()
            fio_output.close()

            fio_output_writer = open(os.path.join(output_dir, os.path.basename(evadd), 'outfiles',
                                                  'ffproc.ampstt.' + str_i), 'w')

            assert(len(passed_staev_band) == len(fi_output)), 'ERROR in length'

            fio_output_writer.writelines(first_line)
            fio_output_writer.writelines(second_line)

            for j in range(len(fi_output)):
                fi_line = fi_output[j].split(',')
                #del_blnk = float(fi_line[0])
                #sdep_blnk = float(fi_line[1])
                #slat_blnk = float(fi_line[2])
                #az_blnk = float(fi_line[3])
                #tt_blnk = float(fi_line[4])
                #cor_blnk = float(fi_line[6])
                tt_tcor_blnk = float(fi_line[5])

                indx = passed_staev_band[j][0]
                grp = passed_staev_band[j][1]
                xcorr = passed_staev_band[j][2]
                da = passed_staev_band[j][3]
                tobs = passed_staev_band[j][4]
                dt = tt_tcor_blnk
                stlat = passed_staev_band[j][6]
                stlon = passed_staev_band[j][7]
                epi = passed_staev_band[j][8]
                sta_id = passed_staev_band[j][9]
                fio_output_writer.writelines('%s  %s  %s  %s  %s  %s  %s  %s  %s  %s\n'
                                             % (indx, grp, xcorr, da, tobs, dt, stlat, stlon, epi, sta_id))
            print '--------------'

        except Exception, e:
            print 'ERROR: %s' % e
            failed += 1
    return passed_staev



########################################################################
############################# Main Program #############################
########################################################################
bands = sys.argv[1]
bands = range(int(bands[0]), int(bands[-1])+1)
band_period = {'1': 30.0, '2': 21.2, '3': 15.0, '4': 10.6, '5': 7.5, '6': 5.3, '7': 3.7, '8': 2.7}

if not single_event:
    proc_ev_ls = glob.glob(os.path.join(remote_dir, '*.*.*.*'))
    for i in range(len(proc_ev_ls)):
        reader_writer(proc_ev_ls[i], bands, phase)
#else:
#    evadd = sys.argv[2]
#    evname = os.path.baseanme(evadd)
#    if not evname:
#        evname = evadd.split('/')[-2]
#    proc_ev_ls = glob.glob(os.path.join(remote_dir, '*.*.*.*'))
