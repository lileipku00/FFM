#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-----------------------------------------------------------------------
#   Filename:  util_PyYSPEC.py
#   Purpose:   PyYSPEC modules
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import os
import subprocess


#--------------------- create_source_inp -------------------------------
def create_source_inp(indir):
    """
    Create source and yspec.in files
    """
    fio_source = open(os.path.join(indir, 'in.yspec.source'))
    fi_source = fio_source.readlines()
    fio_source.close()
    evla, evlo, evdp = fi_source[0].split('\n')[0].split()
    mrr, mtt, mpp, mrt, mrp, mtp = fi_source[1].split('\n')[0].split()
    evyear, evmon, evday, evhour, evmin, evsec = fi_source[3].split('\n')[0].split()

    fio_inv = open(os.path.join(indir, 'SeismogramInventory'))
    fi_inv = fio_inv.readlines()
    fio_inv.close()

    sta_info = []
    for _i in xrange(2, len(fi_inv)):
        sta_info.append(fi_inv[_i].split(';')[4:12])

    sta_info.sort()
    # XXX not sure about the functionality of del_index!
    del_index = []
    fio_sta_yspec = open(os.path.join(indir, 'sta_yspec'), 'w')
    for _i in xrange(len(sta_info)):
        if not sta_info[_i] == []:
            msg = ','.join(sta_info[_i][0:3]) + ',x00,' + ','.join(sta_info[_i][3:9]) + ',' + \
                  ','.join([evla, evlo, evdp]) + ',\n'
            fio_sta_yspec.writelines(msg)
        else:
            del_index.append(_i)
    fio_sta_yspec.close()

    if del_index:
        for _i in del_index.reverse():
            sta_info.remove(sta_info[_i])

    fio_yspecin = open(os.path.join(indir, 'yspec.in_PyYSPEC'), 'r')
    fi_yspecin = fio_yspecin.readlines()
    fio_yspecin.close()

    search = '# source depth (km)'
    for _i in xrange(len(fi_yspecin)):
        if fi_yspecin[_i].find(search) != -1:
            fi_yspecin[_i+1] = '  ' + evdp + '\n'
            break
    search = '# source latitude (deg)'
    for _i in xrange(len(fi_yspecin)):
        if fi_yspecin[_i].find(search) != -1:
            fi_yspecin[_i+1] = '  ' + evla + '\n'
            break
    search = '# source longitude (deg)'
    for _i in xrange(len(fi_yspecin)):
        if fi_yspecin[_i].find(search) != -1:
            fi_yspecin[_i+1] = '  ' + evlo + '\n'
            break
    search = '# M_{r,r} (Nm)'
    for _i in xrange(len(fi_yspecin)):
        if fi_yspecin[_i].find(search) != -1:
            fi_yspecin[_i+1] = '  ' + mrr + '\n'
            break
    search = '# M_{r,theta} (Nm)'
    for _i in xrange(len(fi_yspecin)):
        if fi_yspecin[_i].find(search) != -1:
            fi_yspecin[_i+1] = '  ' + mrt + '\n'
            break
    search = '# M_{r,phi} (Nm)'
    for _i in xrange(len(fi_yspecin)):
        if fi_yspecin[_i].find(search) != -1:
            fi_yspecin[_i+1] = '  ' + mrp + '\n'
            break
    search = '# M_{theta,theta} (Nm)'
    for _i in xrange(len(fi_yspecin)):
        if fi_yspecin[_i].find(search) != -1:
            fi_yspecin[_i+1] = '  ' + mtt + '\n'
            break
    search = '# M_{theta,phi} (Nm)'
    for _i in xrange(len(fi_yspecin)):
        if fi_yspecin[_i].find(search) != -1:
            fi_yspecin[_i+1] = '  ' + mtp + '\n'
            break
    search = '# M_{phi,phi} (Nm)'
    for _i in xrange(len(fi_yspecin)):
        if fi_yspecin[_i].find(search) != -1:
            fi_yspecin[_i+1] = '  ' + mpp + '\n'
            break
    search = '# number of receivers'
    for _i in xrange(len(fi_yspecin)):
        if fi_yspecin[_i].find(search) != -1:
            fi_yspecin[_i+1] = '  ' + str(len(sta_info)) + '\n'
            break

    receivers = []
    for _i in xrange(len(sta_info)):
        receivers.append('  ' + sta_info[_i][4] + '  ' + sta_info[_i][5] + '\n')

    search = '# receiver latitudes and longitudes'
    for _i in xrange(len(fi_yspecin)):
        if fi_yspecin[_i].find(search) != -1:
            fi_yspecin[_i+2:] = receivers
            break

    fio_yspecin_new = open(os.path.join(indir, 'yspec.in'), 'w')
    for _i in xrange(len(fi_yspecin)):
        fio_yspecin_new.write(fi_yspecin[_i])
    fio_yspecin_new.close()


#--------------------- run_yspec ------------------------------------------
def run_yspec(indir_submit, yspec_inp, num_proc, indir_output):
    """
    run yspec
    """
    current_dir = os.getcwd()
    os.chdir(os.path.join(indir_submit))
    subprocess.check_call(['./submit.sh', yspec_inp, num_proc, '-o', indir_output])
    os.chdir(current_dir)
