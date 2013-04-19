import glob
import os


pdata_add = '/import/neptun-radler/AmplitudeProjects/pdata_processed/psdata_events'

evs_ls = glob.glob(os.path.join(pdata_add, '*.*.*.*'))
print '%s events found in the archive!' %(len(evs_ls))
print '\nProblematic events:'
first_line = '#eventID,lat,lon,depth,mag,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp\n'
ev_info = []
for i in range(len(evs_ls)):
    try:
        ev = evs_ls[i]
        ev_name = ev.split('/')[-1]
        fio_source = open(os.path.join(ev, 'outfiles', 'ampinv.source'), 'r')
        f_source = fio_source.readlines()
        ev_year, ev_julianday, ev_hr, ev_min, ev_sec, ev_msec = f_source[1].split()
        evlat, evlon, catalog_depth, inverted_depth = f_source[3].split()
        try:
            mrr, mtt, mpp, mrt, mrp, mtp = f_source[13].split()
        except Exception, e:
            mrr, mtt, mpp, mrt, mrp, mtp = f_source[7].split()
        fio_mag = open(os.path.join(ev, 'README'), 'r')
        f_mag = fio_mag.readlines()
        ev_mag = f_mag[1].split()[1]
        ev_info.append('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n' %(ev_name,evlat,
            evlon,inverted_depth,ev_mag,mrr,mtt,mpp,mrt,mrp,mtp))
    except Exception, e:
        print ev_name

ev_info.insert(0, first_line) 
fio_ls_event = open(os.path.join('.', 'results', 'pdata_events.txt'), 'w')
for i in ev_info:
    fio_ls_event.writelines(i)
fio_ls_event.close() 
