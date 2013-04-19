import os
import shutil
import subprocess

req_phase = 'Pdiff'
req_solver = 'yspec'
req_processes = '32'

add_event_info = '/import/neptun-radler/hosseini-downloads/KASRA/SCRIPTS/gitHUB/myrepo_gitHUB/FFM/results'
add_proc_ev = '/import/neptun-radler/hosseini-downloads/KASRA/FFM'
add_runff = '/home/hosseini/FFINVERSION/AMPLITUDES/Programs/ffproc/source_code'

fio_event_info = open(os.path.join(add_event_info, 'selected_events.txt'))
event_info = fio_event_info.readlines()

#for i in range(len(event_info)):
for i in range(0, 2):
    event_name = event_info[i].split(',')[0]
    event_process = os.path.join(add_proc_ev, event_name)
    if os.path.isdir(event_process):
        shutil.rmtree(event_process)
    os.mkdir(event_process)
    os.chdir(event_process)
    subprocess.check_call(['/home/hosseini/FFINVERSION/AMPLITUDES/Programs/ffproc/FFsetup', req_phase, req_solver, req_processes])
    subprocess.check_call(['matlab', '-nodisplay', '-nosplash', '-nodesktop', '-r', "run('RunFFProcessing')"])
