#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  run_serial.py
#   Purpose:   run the selected events in a serial way
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import os
import shutil
import subprocess
import sys

# ------------------- INPUT -----------------------------
req_phase = 'PP'
req_solver = 'yspec'
req_processes = '40'
update_all = False
add_event_info = '/import/neptun-radler/hosseini-downloads/KASRA/SCRIPTS/gitHUB/myrepo_gitHUB/FFM/INITIALIZATION/results'
add_proc_ev = '/import/neptun-helles/hosseini/FFM_RESULTS/PP_measure_1_sec_LAMBDA_1-5_60_180'
add_runff = '/home/hosseini/FFINVERSION/AMPLITUDES/Programs/ffproc/FFsetup'
# -------------------------------------------------------

exit_flag = raw_input('Did you set run_serial=1 in RunFFProcessing? (y/n)\n')
if exit_flag.upper() != 'Y':
    sys.exit('Please do it first! it is at the top of the file in input section!')

fio_event_info = open(os.path.join(add_event_info, 'selected_events_ALL.txt'))
event_info = fio_event_info.readlines()

for i in range(len(event_info)):
    event_name = event_info[i].split(',')[0]
    event_process = os.path.join(add_proc_ev, event_name)
    if os.path.isdir(event_process):
        if update_all:
            shutil.rmtree(event_process)
        else:
            print "%s has been measured before..." % event_info[i]
            continue
    os.mkdir(event_process)
    os.chdir(event_process)
    subprocess.check_call([add_runff, req_phase, req_solver, req_processes])
    subprocess.check_call(['matlab', '-nodisplay', '-nosplash', '-nodesktop', '-r', "run('RunFFProcessing')"])
