import glob
import os

syn_arch_add = '/import/neptun-helles/hosseini/FFM/YSPEC_SYN_GALLERY'
syn_arch = glob.glob(os.path.join(syn_arch_add, '*.*.*.*'))
for i in range(len(syn_arch)):
    syn_arch[i] = syn_arch[i].split('/')[-1]

selected_event = '/import/neptun-radler/hosseini-downloads/KASRA/SCRIPTS/gitHUB/myrepo_gitHUB/FFM/results/selected_events_all.txt'
updated_event = '/import/neptun-radler/hosseini-downloads/KASRA/SCRIPTS/gitHUB/myrepo_gitHUB/FFM/results/selected_events_updated.txt'
fio_selected_event = open(selected_event, 'r')
fi_selected_event = fio_selected_event.readlines()
fio_selected_event.close()

processed_address = '/import/neptun-helles/hosseini/FFM'
processed_event = glob.glob(os.path.join(processed_address, '*.*.*.*'))

for i in range(len(fi_selected_event)-1, -1, -1):
    try:
        flag_exist = False
        flag_processed = False
        for j in range(len(syn_arch)):
            if fi_selected_event[i].split(',')[0] == syn_arch[j]:
                del fi_selected_event[i]
                break
    except Exception, e:
        print e

fio_selected_event = open(updated_event, 'w')
for i in range(len(fi_selected_event)):
    fio_selected_event.writelines(fi_selected_event[i])
fio_selected_event.close()
