"""
Function to calculate the backazimuth and R, T components of horizontal
components
"""
import glob
import os

from obspy import read
from obspy.core.util import gps2DistAzimuth
from obspy.signal import rotate


# ---------------- INPUT
list_events = '/import/neptun-radler/hosseini-downloads/KASRA/SCRIPTS/' \
              'gitHUB/myrepo_gitHUB/FFM/INITIALIZATION/results/' \
              'selected_events_ALL.txt'
add_syn_data = '/import/neptun-helles/hosseini/SYN_GALLERIES/YSPEC_SYN_GALLERY_1hz'
selected_events = ['0274.2009.273.a']
# ---------------- END INPUT

# Read select_events_ALL.txt file
fio_ls_evs = open(list_events, 'r')
ls_evs = fio_ls_evs.readlines()

for ev in ls_evs:
    try:
        if not selected_events == 'all':
            if not ev.split(',')[0] in selected_events:
                continue
        add_ev = os.path.join(add_syn_data, ev.split(',')[0])
        ls_E = glob.glob(os.path.join(add_ev, 'SAC_realName', '*E'))
        for i in range(len(ls_E)):
            try:
                st_E = read(ls_E[i], format='SAC')[0]
                st_N = read(os.path.join(add_ev, 'SAC_realName',
                                         'grf.%s.%s.%s.x00.%sN'
                                         % (st_E.stats.network,
                                            st_E.stats.station,
                                            st_E.stats.location,
                                            st_E.stats.channel[:-1])),
                            format='SAC')[0]

                (dist, azi, bazi) = gps2DistAzimuth(st_E.stats.sac.evla,
                                                    st_E.stats.sac.evlo,
                                                    st_E.stats.sac.stla,
                                                    st_E.stats.sac.stlo)

                (tr_data_R, tr_data_T) = rotate.rotate_NE_RT(st_N.data,
                                                             st_E.data,
                                                             bazi)

                tr_R = st_N.copy()
                tr_T = st_N.copy()

                tr_R.data = tr_data_R
                tr_R.stats.channel = 'BHR'
                tr_T.data = tr_data_T
                tr_T.stats.channel = 'BHT'
                #tr_R.write(os.path.join(add_ev, 'BH',
                #                        'dis.%s.%s.%s'
                #                        % (tr_R.stats.station,
                #                           tr_R.stats.location,
                #                           tr_R.stats.channel)),
                #           format='SAC')
                tr_T.write(os.path.join(add_ev, 'SAC_realName',
                                         'grf.%s.%s.%s.x00.%s'
                                         % (tr_T.stats.network,
                                            tr_T.stats.station,
                                            tr_T.stats.location,
                                            tr_T.stats.channel)),
                           format='SAC')
            except Exception, e:
                print '%s: %s' % (ls_E[i], e)
    except Exception, e:
        print e
