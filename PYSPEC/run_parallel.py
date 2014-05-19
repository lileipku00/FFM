#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  run_parallel.py
#   Purpose:   run the selected events in parallel
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import multiprocessing
import os
import shutil
import subprocess
import sys
import time

# ------------------- INPUT -----------------------------
req_phase = 'Pdiff'
req_solver = 'yspec'
req_processes = '1'
add_event_info = '/import/neptun-radler/hosseini-downloads/KASRA/SCRIPTS/gitHUB/myrepo_gitHUB/FFM/INITIALIZATION/results'
add_proc_ev = '/import/neptun-helles/hosseini/FFM_RESULTS/Pdiff_measure_1_sec_LAMBDA_1-5_90_180_TEST'
add_runff = '/home/hosseini/FFINVERSION/AMPLITUDES/Programs/ffproc/FFsetup'

# Number of parallel requests (PLEASE DO NOT USE MORE THAN 5...License issue!)
req_np = 5
# -------------------------------------------------------


##################### FUNCTIONS #########################
def serial_run(add_proc_ev, event_info, add_runff, req_phase, req_solver, req_processes):
    """
    One run: FFsetup + RunFFprocessing
    This will be distributed over the number of required processes
    """
    event_name = event_info.split(',')[0]
    event_process = os.path.join(add_proc_ev, event_name)
    if os.path.isdir(event_process):
        shutil.rmtree(event_process)
    os.mkdir(event_process)
    os.chdir(event_process)
    subprocess.check_call([add_runff, req_phase, req_solver, req_processes])
    subprocess.check_call(['matlab', '-nodisplay', '-nosplash', '-nodesktop', '-r', "run('RunFFProcessing')"])

##################### MAIN PROGRAM #########################
exit_flag = raw_input('Did you set run_serial=1 in RunFFProcessing? (y/n)\n')
if exit_flag.upper() != 'Y':
    sys.exit('Please do it first! it is at the top of the file in input section!')

fio_event_info = open(os.path.join(add_event_info, 'selected_events.txt'))
event_info = fio_event_info.readlines()

parallel_len_run = range(0, len(event_info))
len_par_grp = [parallel_len_run[n:n+req_np] for n in range(0, len(parallel_len_run), req_np)]
par_jobs = []
for i in range(len(event_info)):
    p = multiprocessing.Process(target=serial_run, args=(add_proc_ev, event_info[i], add_runff, req_phase, req_solver,
                                                         req_processes))
    par_jobs.append(p)
for l in range(len(len_par_grp)):
    for ll in len_par_grp[l]:
        par_jobs[ll].start()
        time.sleep(5)
    par_jobs[ll].join()
