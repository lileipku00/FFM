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
req_phase = 'Pdiff'
req_solver = 'yspec'
req_processes = '64'
add_event_info = '/import/neptun-radler/hosseini-downloads/KASRA/SCRIPTS/gitHUB/myrepo_gitHUB/FFM/results'
add_proc_ev = '/import/neptun-radler/hosseini-downloads/KASRA/FFM'
add_runff = '/home/hosseini/FFINVERSION/AMPLITUDES/Programs/ffproc/FFsetup'
# -------------------------------------------------------

exit_flag = raw_input('Did you uncomment the Exit option in RunFFProcessing? (y/n)\n')
if exit_flag.upper() != 'Y':
    sys.exit('Please uncomment the flag at the very end of the file!')

fio_event_info = open(os.path.join(add_event_info, 'selected_events.txt'))
event_info = fio_event_info.readlines()

for i in range(len(event_info)):
    event_name = event_info[i].split(',')[0]
    event_process = os.path.join(add_proc_ev, event_name)
    if os.path.isdir(event_process):
        shutil.rmtree(event_process)
    os.mkdir(event_process)
    os.chdir(event_process)
    subprocess.check_call([add_runff, req_phase, req_solver, req_processes])
    subprocess.check_call(['matlab', '-nodisplay', '-nosplash', '-nodesktop', '-r', "run('RunFFProcessing')"])
