import glob
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime
import os
import sys

####################### output_file_reader #############################


def output_file_reader(evinfo, req_band='band01'):
    """
    This function reads one output file and pass all the values for further analysis in next steps
    """
    add_output = evinfo[0]
    empty_array = np.array([])
    if not os.path.isfile(os.path.join(add_output, 'outfiles', 'ffproc.ampstt.%s' % req_band)):
        print "%s is not found!" % os.path.join(add_output, 'outfiles', 'ffproc.ampstt.%s' % req_band)
        return empty_array
    if not os.path.isfile(os.path.join(add_output, 'outfiles', 'ffproc.receivers')):
        print "%s is not found!" % os.path.join(add_output, 'outfiles', 'ffproc.receivers')
        return empty_array
    rd_output = np.loadtxt(os.path.join(add_output, 'outfiles', 'ffproc.ampstt.%s' % req_band), dtype='S', comments='#')
    sta_output = np.loadtxt(os.path.join(add_output, 'outfiles', 'ffproc.receivers'), dtype='S', comments='#')
    new_col_cr = np.empty([np.shape(rd_output)[0], 9], dtype=object)
    new_col_cr[:, 0] = (sta_output[:, 5].astype(np.float) - sta_output[:, 6].astype(np.float))/1000.
    new_col_cr[:, 1] = os.path.basename(add_output)
    new_col_cr[:, 2] = evinfo[1]
    new_col_cr[:, 3] = evinfo[2]
    new_col_cr[:, 4] = evinfo[3]
    new_col_cr[:, 5] = evinfo[4]
    new_col_cr[:, 6] = evinfo[5]
    new_col_cr[:, 7] = evinfo[6]
    new_col_cr[:, 8] = req_band
    output_sta_evname = np.append(rd_output, new_col_cr, 1)
    return output_sta_evname

####################### event_filter #############################


def event_filter(par_add, selected_events_add, all_events=True, min_dp=-10, max_dp=1000):
    """
    Filters the events in one par_dir based on the required inputs
    """
    if not os.path.isdir(par_add):
        print "%s is not a valid directory!" % par_add
        return False, False

    selected_events = np.loadtxt(fname=selected_events_add, dtype='S', comments='#', delimiter=',')
    event_adds = []
    for i in range(len(selected_events[:, 1])):
        if all_events != True:
            if not selected_events[:, 1][i] in all_events: continue
        event_adds.append([os.path.join(par_add, selected_events[:, 1][i]), selected_events[:, 0][i]])
    passed_event_adds = []
    for i in range(len(event_adds)):
        add_flag = False
        try:
            fio_source = open(os.path.join(event_adds[i][0], 'outfiles', 'ffproc.source'), 'r')
            f_source = fio_source.readlines()
            ev_year, ev_julianday, ev_hr, ev_min, ev_sec, ev_msec = f_source[1].split()
            evlat, evlon, catalog_depth, inverted_depth = f_source[3].split()
            try:
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[13].split()
            except Exception, e:
                print 'WARNING: inverted moment tensor was not found: %s' % event_adds[i][0]
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[7].split()
        except Exception, e:
            print 'ERROR: %s' % e
        ev_date = UTCDateTime(year=int(ev_year), julday=int(ev_julianday))
        ev_date_str = '%4s%2s%2s' % (ev_date.year, ev_date.month, ev_date.day)
        ev_date_str = ev_date_str.replace(' ', '0')
        ev_time = '%2s%2s%2s' % (ev_hr, ev_min, ev_sec)
        ev_time = ev_time.replace(' ', '0')
        ev_id = event_adds[i][1]

        # Check for depth
        if min_dp <= float(inverted_depth) < max_dp:
            add_flag = True
        if add_flag:
            passed_event_adds.append([event_adds[i][0], ev_date_str, ev_time, ev_id, evlat, evlon, inverted_depth])
    return passed_event_adds

####################### array_station_filter #############################


def array_station_filter(passed_array, min_xcorr=-100, max_xcorr=100, min_epi=0., max_epi=360., check_clip=True):
    """
    Filters the stations in an array based on the required inputs
    """
    empty_array = np.array([])
    pass_stas_corr_1 = passed_array[passed_array[:, 6].astype(np.float) >= min_xcorr]

    if not pass_stas_corr_1.size == 0:
         pass_stas_corr_2 = pass_stas_corr_1[pass_stas_corr_1[:, 6].astype(np.float) < max_xcorr]
    else:
        return empty_array

    passed_stas_epi_1 = pass_stas_corr_2[pass_stas_corr_2[:, 4].astype(np.float) >= min_epi]
    if not passed_stas_epi_1.size == 0:
        passed_stas_epi_2 = passed_stas_epi_1[passed_stas_epi_1[:, 4].astype(np.float) < max_epi]
    else:
        return empty_array
    if passed_stas_epi_2.size == 0:
        return empty_array

    if check_clip:
        passed_stas = passed_stas_epi_2[passed_stas_epi_2[:, 19].astype(np.float) < 0.1]

    return passed_stas

####################### check_selection #############################


def check_selection(filt_array, min_xcorr, max_xcorr, min_epi, max_epi, check_clip=True):
    """
    Check whether the selection procedure worked well
    """
    plt.subplot(2, 2, 1)
    plt.plot(filt_array[:, 2].astype(np.float), filt_array[:, 6].astype(np.float), 'r.', linewidth=15)
    plt.xlabel('Latitude', size='large', weight='bold')
    plt.ylabel('xcorr factor', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.hlines(min_xcorr, np.min(filt_array[:, 2].astype(np.float)), np.max(filt_array[:, 2].astype(np.float)), 'k',
               linestyles='--')
    plt.hlines(max_xcorr, np.min(filt_array[:, 2].astype(np.float)), np.max(filt_array[:, 2].astype(np.float)), 'k',
               linestyles='--')
    plt.subplot(2, 2, 2)
    plt.plot(filt_array[:, 3].astype(np.float), filt_array[:, 6].astype(np.float), 'r.', linewidth=15)
    plt.xlabel('Longitude', size='large', weight='bold')
    plt.ylabel('xcorr factor', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.hlines(min_xcorr, np.min(filt_array[:, 3].astype(np.float)), np.max(filt_array[:, 3].astype(np.float)), 'k',
               linestyles='--')
    plt.hlines(max_xcorr, np.min(filt_array[:, 3].astype(np.float)), np.max(filt_array[:, 3].astype(np.float)), 'k',
               linestyles='--')
    plt.subplot(2, 2, 3)
    plt.plot(filt_array[:, 4].astype(np.float), filt_array[:, 6].astype(np.float), 'b.', linewidth=15)
    plt.xlabel('Epicentral Distance', size='large', weight='bold')
    plt.ylabel('xcorr factor', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.vlines(min_epi, min_xcorr, max_xcorr, 'k', linestyles='--')
    plt.vlines(max_epi, min_xcorr, max_xcorr, 'k', linestyles='--')
    plt.subplot(2, 2, 4)
    plt.plot(filt_array[:, 19].astype(np.float), filt_array[:, 6].astype(np.float), 'r.', linewidth=15)
    plt.xlabel('Clip value', size='large', weight='bold')
    plt.ylabel('xcorr factor', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.show()

####################### raydata_input_generator #############################


def raydata_input_generator(filt_array, input_file_name, twinned, phase, min_xcorr, min_depth, max_depth,
                            min_epi, max_epi, check_clip):
    """
    Generate input file compatible with raydata input files
    """
    if not os.path.isfile(os.path.join(os.path.curdir, 'info', 'Pdef_%s' % phase)):
        sys.exit('%s could not be found!' % os.path.join(os.path.curdir, 'info', 'Pdef_%s' % phase))
    phase_def = open(os.path.join(os.path.curdir, 'info', 'Pdef_%s' % phase), 'r').readlines()
    if not os.path.isfile(os.path.join(os.path.curdir, 'info', 'bpf.omega_m')):
        sys.exit('%s could not be found!' % os.path.join(os.path.curdir, 'info', 'bpf.omega_m'))
    filt_def = open(os.path.join(os.path.curdir, 'info', 'bpf.omega_m'), 'r').readlines()

    inp_lines = []
    inp_lines.append('%s\n' % input_file_name)
    inp_lines.append('%s\n' % twinned)
    inp_lines.append('# Phase : %s\n' % phase)
    inp_lines.append('# Traveltime data, xcorr_min : %s\n' % min_xcorr)
    inp_lines.append('# depth_min : %s, depth_max : %s\n' % (min_depth, max_depth))
    inp_lines.append('# epi_min : %s, epi_max : %s\n' % (min_epi, max_epi))
    inp_lines.append('# Clip : %s\n' % check_clip)
    for ln in phase_def:
        inp_lines.append(ln)
    for ln in filt_def:
        inp_lines.append(ln)
    for sta in filt_array:
        netid = sta[5].split('.')[0]
        staid = sta[5].split('.')[1]
        chaid = sta[5].split('.')[4]
        first_line = '%s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s 1  0  0\n' % (sta[24], sta[25], sta[26],
                                                                                     sta[1], staid, netid, chaid,
                                                                                     sta[27], sta[28], sta[29],
                                                                                     sta[2], sta[3], sta[22])
        second_line = '1    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0\n'
        third_line = '      %s  %s  %s  1  %s\n' % (sta[7], sta[9], sta[6], sta[18])
        inp_lines.append(first_line)
        inp_lines.append(second_line)
        inp_lines.append(third_line)
    input_file_fio = open(input_file_name, 'w')
    input_file_fio.writelines(inp_lines)

####################### raydata_input #############################


def raydata_input(bg_model, input_file_name, phase):
    """
    make in.input_file_name for raydata
    """
    in_input_fio = open('in.%s' % input_file_name, 'w')
    in_input_fio.write('%s\n' % bg_model)
    in_input_fio.write('1  20\n')
    in_input_fio.write('%s\n' % input_file_name)
    in_input_fio.write('Pdef_%s' % phase)
    in_input_fio.close()



####################### station_filter #############################


def station_filter(ls_stas, min_xcorr=-100, max_xcorr=100, min_epi=0., max_epi=360., check_clip=True):
    """
    Filters the stations based on the required inputs
    """
    empty_array = np.array([])
    pass_stas_corr_1 = ls_stas[min_xcorr <= ls_stas[:, 6].astype(np.float)]

    if not pass_stas_corr_1.size == 0:
         pass_stas_corr_2 = pass_stas_corr_1[max_xcorr > pass_stas_corr_1[:, 6].astype(np.float)]
    else:
        return empty_array

    passed_stas_epi_1 = pass_stas_corr_2[min_epi <= pass_stas_corr_2[:, 4].astype(np.float)]
    if not passed_stas_epi_1.size == 0:
        passed_stas_epi_2 = passed_stas_epi_1[max_epi > passed_stas_epi_1[:, 4].astype(np.float)]
    else:
        return empty_array
    if passed_stas_epi_2.size == 0:
        return empty_array

    if check_clip:
        passed_stas = passed_stas_epi_2[passed_stas_epi_2[:, 19].astype(np.float) < 0.1]

    return passed_stas