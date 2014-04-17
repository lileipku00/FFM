import glob
import os

syn_arch_add = '/import/neptun-helles/hosseini/FFM/YSPEC_SYN_GALLERY'
syn_arch = glob.glob(os.path.join(syn_arch_add, '*.*.*.*'))
for i in range(len(syn_arch)):
    syn_arch[i] = syn_arch[i].split('/')[-1]

selected_event = '/import/neptun-radler/hosseini-downloads/KASRA/SCRIPTS/gitHUB/myrepo_gitHUB/FFM/INITIALIZATION/results/selected_events.txt'
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
                flag_exist = True
        if flag_exist:
            for j in range(len(processed_event)):
                if fi_selected_event[i].split(',')[0] == processed_event[j].split('/')[-1]:
                    flag_processed = True
        if flag_exist:
            if flag_processed:
                del fi_selected_event[i]
            else:
                continue
        else:
            del fi_selected_event[i]
    except Exception, e:
        print e

fio_selected_event = open(selected_event, 'w')
for i in range(len(fi_selected_event)):
    fio_selected_event.writelines(fi_selected_event[i])
fio_selected_event.close()
