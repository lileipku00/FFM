import glob
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime
import os
import shutil
import sys
import time
import pyvtk as pvtk

####################### output_file_reader #############################


def output_file_reader(evinfo, req_band='band01'):
    """
    This function reads one output file and pass all the values for further analysis in next steps
    ATTENTION: there is not filtering at this stage, it reads all...
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
    # Generate an empty array with size of: rd_output_rows X 9
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
    # Now append the rd_output_rows X 9 array to the original rd_output
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
        except Exception, e:
            print 'ERROR: %s' % e

    return passed_event_adds

####################### array_station_filter #############################


def array_station_filter(passed_array, min_xcorr=-100, max_xcorr=100, min_epi=0., max_epi=360., check_clip=True):
    """
    Filters the stations in an array based on the required inputs
    """
    empty_array = np.array([])

    # --------------- XCORRELATION ---------------
    pass_stas_corr_1 = passed_array[passed_array[:, 6].astype(np.float) >= min_xcorr]
    if not pass_stas_corr_1.size == 0:
        pass_stas_final = pass_stas_corr_1[pass_stas_corr_1[:, 6].astype(np.float) < max_xcorr]
    else:
        return empty_array

    # --------------- EPICENTRAL ---------------
    passed_stas_epi_1 = pass_stas_final[pass_stas_final[:, 4].astype(np.float) >= min_epi]
    if not passed_stas_epi_1.size == 0:
        passed_stas_final = passed_stas_epi_1[passed_stas_epi_1[:, 4].astype(np.float) < max_epi]
    else:
        return empty_array

    # --------------- CHECK CLIPS ---------------
    if check_clip:
        passed_stas_final = passed_stas_final[passed_stas_final[:, 19].astype(np.float) < 0.1]

    if not passed_stas_final.size == 0:
        return passed_stas_final
    else:
        return empty_array

####################### array_station_filter_mark #############################


def array_station_filter_mark(all_output_files, min_xcorr=-100, max_xcorr=100, min_epi=0., max_epi=360.,
                              check_clip=True):
    """
    Filters the stations in an array based on the required inputs
    """
    # --------------- XCORRELATION ---------------
    all_output_files[:, 10, :][all_output_files[:, 6, :].astype(np.float) < min_xcorr] = -1
    all_output_files[:, 10, :][all_output_files[:, 6, :].astype(np.float) >= max_xcorr] = -1

    # --------------- EPICENTRAL ---------------
    all_output_files[:, 10, :][all_output_files[:, 4, :].astype(np.float) < min_epi] = -1
    all_output_files[:, 10, :][all_output_files[:, 4, :].astype(np.float) >= max_epi] = -1

    # --------------- CHECK CLIPS ---------------
    if check_clip:
        all_output_files[:, 10, :][all_output_files[:, 19, :].astype(np.float) > 0.1] = -1

    return all_output_files

####################### check_selection #############################


def check_selection(filt_array, min_xcorr, max_xcorr, min_epi, max_epi, check_clip=True):
    """
    Check whether the selection procedure worked well
    """
    plt.ion()
    plt.figure()
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

####################### axi_kernel_receiver_writer #############################


def axi_kernel_receiver_writer(evstas):
    """
    Create AXISEM receiver.dat file for further analysis using AXISEM
    """
    print '\n======>> Create AXISEM receiver.dat file for Kernel calculations'
    if not os.path.isdir(os.path.join(os.path.curdir, 'RESULTS', evstas[0, 23, 0])):
        os.mkdir(os.path.join(os.path.curdir, 'RESULTS', evstas[0, 23, 0]))
    else:
        print 'Directory already exists: %s' % os.path.join(os.path.curdir, 'RESULTS', evstas[0, 23, 0])
        return

    band_period = {'band01': 30.0, 'band02': 21.2, 'band03': 15.0, 'band04': 10.6, 'band05': 7.5, 'band06': 5.3,
                   'band07': 3.7, 'band08': 2.7}

    line_write = []
    counter = 0
    for j in range(evstas.shape[0]):
        len_freq_avail = len(evstas[j, :, evstas[j, 10, :] != -1])
        if len_freq_avail == 0:
            continue
        third_line = '%s_%s  %s  %s  %s\n' % (evstas[j, 5, 0].split('.')[0], evstas[j, 5, 0].split('.')[1],
                                              evstas[j, 2, 0], evstas[j, 3, 0], len_freq_avail)
        line_write.append(third_line)
        for i in range(len(evstas[j, 10, :])):
            if evstas[j, 10, i] == -1:
                continue
            period = band_period[evstas[j, 30, i]]
            forth_line = 'kernel_%s  Gabor_%s  CC  %s  %s\n' \
                         % (str(period).replace('.', '_'), period, evstas[j, 17, i],
                            float(evstas[j, 17, i])+float(evstas[j, 18, i]))
            line_write.append(forth_line)
        counter += 1
    first_line = '%i\n' % counter
    second_line = 'vp  Z\n'
    receiver_fio = open(os.path.join(os.path.curdir, 'RESULTS', evstas[0, 23, 0], 'receiver.dat'), 'w')
    receiver_fio.writelines(first_line)
    receiver_fio.writelines(second_line)
    for ln in line_write:
        receiver_fio.writelines(ln)
    receiver_fio.close()

    readme_fio = open(os.path.join(os.path.curdir, 'RESULTS', evstas[0, 23, 0], 'README.txt'), 'w')
    readme_fio.writelines('number_of_receivers: %i\n' % counter)
    readme_fio.close()

    #plt.ion()
    #plt.figure()
    #plt.subplot(2, 1, 1)
    #plt.plot(filt_array[:, 2].astype(np.float), t_corr, 'r.')
    #plt.xlabel('Latitude', size='large', weight='bold')
    #plt.ylabel('Common Correction Values', size='large', weight='bold')
    #plt.xticks(size='large', weight='bold')
    #plt.yticks(size='large', weight='bold')
    #plt.subplot(2, 1, 2)
    #plt.plot(filt_array[:, 3].astype(np.float), t_corr, 'r.')
    #plt.xlabel('Longitude', size='large', weight='bold')
    #plt.ylabel('Common Correction Values', size='large', weight='bold')
    #plt.xticks(size='large', weight='bold')
    #plt.yticks(size='large', weight='bold')
    #plt.show

####################### axi_kernel_receiver_combiner #############################


def axi_kernel_receiver_combiner(measure_par_add, measures):
    """
    For several type of measurements (P, Pdiff, ...), combine the receiver.dat file
    Usage:
    from output_reader import axi_kernel_receiver_combiner
    axi_kernel_receiver_combiner('./RESULTS', ['P', 'Pdiff'])
    """
    ls_all_events = np.array([], dtype='object')
    for i in range(len(measures)):
        ls_all_events = np.append(ls_all_events, glob.glob(os.path.join(measure_par_add, measures[i], '*.*.*.*')))
    for i in range(len(ls_all_events)):
        ls_all_events[i] = os.path.basename(ls_all_events[i])
    ls_all_events_unique = np.unique(ls_all_events)
    for i in range(len(ls_all_events_unique)):
        cont_flag = True
        for j in range(len(measures)):
            if not os.path.isfile(os.path.join(measure_par_add, measures[j], ls_all_events_unique[i], 'receiver.dat')):
                cont_flag = False
        if not cont_flag:
            continue
        all_staevs = []
        len_staev = 0
        for j in range(len(measures)):
            receiver_tmp_fio = open(os.path.join(measure_par_add, measures[j], ls_all_events_unique[i],
                                                 'receiver.dat'), 'r')
            receiver_tmp = receiver_tmp_fio.readlines()
            all_staevs.append(receiver_tmp[2:])

            len_staev_fio = open(os.path.join(measure_par_add, measures[j], ls_all_events_unique[i],
                                              'README.txt'), 'r')
            len_staev += int(len_staev_fio.readlines()[0].split(':')[1])
        receiver_all_fio = open(os.path.join(measure_par_add,
                                             'receiver_%s.dat' % ls_all_events_unique[i].replace('.', '_')),
                                'w')
        receiver_all_fio.writelines('%i\n' % len_staev)
        receiver_all_fio.writelines('vp   Z\n')
        for ln_all in all_staevs:
            receiver_all_fio.writelines(ln_all)
        receiver_all_fio.close()

    print '\n================================='
    print 'WARNING: vp and Z are hard coded!'
    print '================================='





####################### compile_raydata_raymatrix #############################


def compile_raydata_raymatrix():
    """
    Compile both raydata and raymatrix for further usage
    """
    cur_dir = os.path.abspath(os.curdir)
    os.chdir(os.path.join(os.curdir, 'raydata_raymatrix', 'raydata_src'))
    os_sys = os.system('./make')
    if not os_sys == 0:
        print "raydata can not be compiled properly"
    os.chdir(cur_dir)

    os.chdir(os.path.join(os.curdir, 'raydata_raymatrix', 'raymatrix_src'))
    os_sys = os.system('./make')
    if not os_sys == 0:
        print "raydata can not be compiled properly"
    os.chdir(cur_dir)

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
    input_file_fio = open(os.path.join('RESULTS', input_file_name), 'w')
    input_file_fio.writelines(inp_lines)

####################### raydata_input #############################


def raydata_input(bg_model, input_file_name, phase, max_num_arrival, delay_wrt_first_arrival):
    """
    make in.input_file_name for raydata
    """
    in_input_fio = open(os.path.join('RESULTS', 'in.raydata_%s' % input_file_name), 'w')
    in_input_fio.write('%s\n' % bg_model)
    in_input_fio.write('%s  %s\n' % (max_num_arrival, delay_wrt_first_arrival))
    in_input_fio.write('%s\n' % input_file_name)
    in_input_fio.write('Pdef_%s' % phase)
    in_input_fio.close()

####################### raymatrix_input #############################


def raymatrix_input(vp_vs_Qs, kernel_quad_km, vertex_file, facet_file, input_file_name):
    """
    create in.raymatrix_input_file_name
    """
    in_input_fio = open(os.path.join('RESULTS', 'in.raymatrix_%s' % input_file_name), 'w')
    in_input_fio.write('%s %s %s\n' % (vp_vs_Qs[0], vp_vs_Qs[1], vp_vs_Qs[2]))
    in_input_fio.write('%s\n' % kernel_quad_km)
    in_input_fio.write('%s\n' % vertex_file)
    in_input_fio.write('%s\n' % facet_file)
    in_input_fio.write('%s' % input_file_name)
    in_input_fio.close()

    in_input_fio = open(os.path.join('RESULTS', 'in.matrixT.%s' % input_file_name), 'w')
    in_input_fio.write('matrixT.%s' % input_file_name)
    in_input_fio.close()

####################### prepare_dir #############################


def prepare_dir(input_file_name):
    """
    prepare directory for one run of raydata and raymatrix
    """
    cur_dir = os.path.abspath(os.curdir)
    if os.path.isdir(os.path.join(os.path.curdir, 'RESULTS', '%s_dir' % input_file_name)):
        sys.exit("Directory already exists!")
    os.mkdir(os.path.join(os.path.curdir, 'RESULTS', '%s_dir' % input_file_name))
    shutil.copy(os.path.join('RESULTS', input_file_name),
                os.path.join(os.path.curdir, 'RESULTS', '%s_dir' % input_file_name))
    shutil.copy(os.path.join('RESULTS', 'in.raydata_%s' % input_file_name),
                os.path.join(os.path.curdir, 'RESULTS', '%s_dir' % input_file_name))
    shutil.copy(os.path.join('RESULTS', 'in.raymatrix_%s' % input_file_name),
                os.path.join(os.path.curdir, 'RESULTS', '%s_dir' % input_file_name))
    shutil.copy(os.path.join('RESULTS', 'in.matrixT.%s' % input_file_name),
                os.path.join(os.path.curdir, 'RESULTS', '%s_dir' % input_file_name))

    files_glob = glob.glob(os.path.join(os.path.curdir, 'raydata_raymatrix', 'files', '*'))
    for fi in files_glob:
        shutil.copy(fi, os.path.join('RESULTS', '%s_dir' % input_file_name))


    shutil.copy(os.path.join(os.curdir, 'raydata_raymatrix', 'raydata_src', 'raydata'),
                os.path.join(cur_dir, 'RESULTS', '%s_dir' % input_file_name))
    shutil.copy(os.path.join(os.curdir, 'raydata_raymatrix', 'raymatrix_src', 'raymatrix'),
                os.path.join(cur_dir, 'RESULTS', '%s_dir' % input_file_name))
    shutil.copy(os.path.join(os.curdir, 'raydata_raymatrix', 'raymatrix_src', 'mat2asc'),
                os.path.join(cur_dir, 'RESULTS', '%s_dir' % input_file_name))

####################### run_raydata_raymatrix #############################


def run_raydata_raymatrix(input_file_name, raydata=True, raymatrix=True):
    """
    run both raydata and raymatrix in the directory
    """
    cur_dir = os.path.abspath(os.curdir)
    if raydata:
        print '\n======>> run raydata at ./RESULTS/%s_dir' % input_file_name
        os.chdir(os.path.join(cur_dir, 'RESULTS', '%s_dir' % input_file_name))
        os_sys = os.system('./raydata < in.raydata_%s' % input_file_name)
        if not os_sys == 0:
            print 'raydata was not executed correctly!'

    if raymatrix:
        print '\n======>> run raymatrix at ./RESULTS/%s_dir' % input_file_name
        os.chdir(os.path.join(cur_dir, 'RESULTS', '%s_dir' % input_file_name))
        os_sys = os.system('./raymatrix < in.raymatrix_%s' % input_file_name)
        if not os_sys == 0:
            print 'raymatrix was not executed correctly!'

    os.chdir(cur_dir)

####################### parallel_raydata_raymatrix #############################


def parallel_raydata_raymatrix(filt_array, input_file_name, twinned, phase, min_xcorr, min_depth, max_depth,
                               min_epi, max_epi, check_clip, bg_model, vp_vs_Qs, kernel_quad_km, vertex_file,
                               facet_file, max_num_arrival, delay_wrt_first_arrival, run_raydata, run_raymatrix):
    """
    To run raydata and raymatrix in parallel
    """
    raydata_input_generator(filt_array=filt_array, input_file_name=input_file_name, twinned=twinned,
                            phase=phase, min_xcorr=min_xcorr, min_depth=min_depth,
                            max_depth=max_depth, min_epi=min_epi, max_epi=max_epi, check_clip=check_clip)
    raydata_input(bg_model=bg_model, input_file_name=input_file_name, phase=phase,
                  max_num_arrival=max_num_arrival, delay_wrt_first_arrival=delay_wrt_first_arrival)
    raymatrix_input(vp_vs_Qs=vp_vs_Qs, kernel_quad_km=kernel_quad_km, vertex_file=vertex_file,
                    facet_file=facet_file, input_file_name=input_file_name)
    print '\n======>> prepare output directory at: ./RESULTS/%s_dir' % input_file_name
    prepare_dir(input_file_name=input_file_name)
    run_raydata_raymatrix(input_file_name=input_file_name, raydata=run_raydata, raymatrix=run_raymatrix)

####################### raydata_ccorr_reader #############################


def raydata_ccorr_reader(filt_array, input_file_name, corr_io_list):
    """
    Read common correction results and rewrite the values in filt_array
    """
    ccorr_arr = np.loadtxt(os.path.join(os.path.curdir, 'RESULTS', '%s_dir' % input_file_name,
                                        'ell_ccor.%s' % input_file_name), dtype='S', comments='#', delimiter=',')
    kd_avail = np.unique(ccorr_arr[:, 0].astype(np.int))
    t_corr = np.zeros(len(ccorr_arr[ccorr_arr[:, 0].astype(np.int) == kd_avail[0]]))
    for i in range(len(kd_avail)):
        tar_ccorr_arr = ccorr_arr[ccorr_arr[:, 0].astype(np.int) == kd_avail[i]]
        if corr_io_list[0] == 1:
            t_corr += tar_ccorr_arr[:, 6].astype(np.float)
        if corr_io_list[1] == 1:
            t_corr += tar_ccorr_arr[:, 7].astype(np.float)
        if corr_io_list[2] == 1:
            t_corr += tar_ccorr_arr[:, 8].astype(np.float)
    plt.ion()
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(filt_array[:, 2].astype(np.float), t_corr, 'r.')
    plt.xlabel('Latitude', size='large', weight='bold')
    plt.ylabel('Common Correction Values', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.subplot(2, 1, 2)
    plt.plot(filt_array[:, 3].astype(np.float), t_corr, 'r.')
    plt.xlabel('Longitude', size='large', weight='bold')
    plt.ylabel('Common Correction Values', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.show()

    filt_array[:, 8] = filt_array[:, 8].astype(np.float) - t_corr
    return filt_array

####################### raydata_ccorr_writer #############################


def raydata_ccorr_writer(filt_array_corr, events_dir):
    """
    Write similar files as ffproc.ampstt.* but with common correction applied
    """
    req_dirs = np.unique(filt_array_corr[:, 23])
    for r_dir in req_dirs:
        os.mkdir(os.path.join(os.path.curdir, 'RESULTS', r_dir))
        os.mkdir(os.path.join(os.path.curdir, 'RESULTS', r_dir, 'outfiles'))
        filt_array_this_dir = filt_array_corr[filt_array_corr[:, 23] == r_dir]
        np.savetxt(os.path.join(os.path.curdir, 'RESULTS', r_dir, 'outfiles',
                                'ffproc.ampstt.%s' % filt_array_this_dir[0, 30]), filt_array_this_dir[:, 0:22],
                   fmt='%s', delimiter='     ')
        shutil.copy(os.path.join(events_dir, r_dir, 'outfiles', 'ffproc.source'),
                    os.path.join(os.curdir, 'RESULTS', r_dir, 'outfiles'))
        shutil.copy(os.path.join(events_dir, r_dir, 'outfiles', 'ffproc.receivers'),
                    os.path.join(os.curdir, 'RESULTS', r_dir, 'outfiles'))
        shutil.copy(os.path.join(events_dir, r_dir, 'outfiles', 'ampinv.source'),
                    os.path.join(os.curdir, 'RESULTS', r_dir, 'outfiles'))

####################### check_par_jobs #############################


def check_par_jobs(jobs, sleep_time=1):
    """
    check whether all the parallel jobs are finished or not
    """
    pp_flag = True
    while pp_flag:
        for proc in jobs:
            if proc.is_alive():
                print '.',
                sys.stdout.flush()
                time.sleep(sleep_time)
                pp_flag = True
                break
            else:
                pp_flag = False
    if not pp_flag:
        print '\n\nAll %s processes are finished...\n' % len(jobs)

####################### mat2asc_run #############################


def mat2asc_run(input_file_name):
    """
    run mat2asc in each directory
    """
    cur_dir = os.path.abspath(os.curdir)
    print '\n======>> run mat2asc at ./RESULTS/%s_dir' % input_file_name
    os.chdir(os.path.join(cur_dir, 'RESULTS', '%s_dir' % input_file_name))
    os_sys = os.system('./mat2asc < in.matrixT.%s' % input_file_name)
    if not os_sys == 0:
        print 'mat2asc was not executed correctly!'
    os.chdir(cur_dir)

####################### vtk_generator #############################


def vtk_generator(input_file_name_part, req_band, vertex_file, facet_file, parallel_exec, len_dirs):
    """
    VTK file generator out of all the results for one complete run
    INFO:
    - facet file is indexed from 0
    - BUT! ascii.* files are indexed from 1, so we need to subtract them (-1)
    - When describing the cells in terms of point indices, the points must be indexed starting at 0.
    """
    print '\n======>> Creating VTK file'
    input_file_name = '%s_%s_%s' % (input_file_name_part, req_band, 1)
    direname = os.path.join(os.path.curdir, 'RESULTS', '%s_dir' % input_file_name)

    print '------> Load Vertex file'
    mesh_points = np.loadtxt(os.path.join(direname, vertex_file), skiprows=2, comments='#')
    print '------> Load Facet file'
    mesh_facets = np.loadtxt(os.path.join(direname, facet_file), dtype=np.int, skiprows=1, comments='#')

    mat_val_all = [0.0]*len(mesh_points)
    for nj in range(len_dirs):
        input_file_name = '%s_%s_%s' % (input_file_name_part, req_band, nj+1)
        direname = os.path.join(os.path.curdir, 'RESULTS', '%s_dir' % input_file_name)
        print '------> create VTK file at ./RESULTS/%s_dir' % input_file_name
        ascii_file = 'ascii.matrixT.%s' % input_file_name
        fmatrix = open(os.path.join(direname, ascii_file), 'r')
        fmatrix_r = fmatrix.readlines()
        mat_indx = []
        mat_val = []
        counter = 0
        for j in range(5, len(fmatrix_r)):
            if counter == 3:
                counter = 0
                continue
            if counter == 0:
                # mat_indx: matrix index
                mat_indx_tmp = fmatrix_r[j].split()
                for i in range(len(mat_indx_tmp)):
                    # IF indexing in ASCII file starts from 0,
                    # we do not need -1.
                    # Otherwise we need it!
                    mat_indx.append(int(mat_indx_tmp[i])-1)
            if counter == 1:
                # mat_val: matrix value
                mat_val_tmp = fmatrix_r[j].split()
                for i in range(len(mat_val_tmp)):
                    mat_val.append(abs(float(mat_val_tmp[i])))
            counter += 1

        for i in range(len(mat_indx)):
            mat_val_all[mat_indx[i]] += mat_val[i]

    vtk = pvtk.VtkData(pvtk.UnstructuredGrid(mesh_points, tetra=mesh_facets),
                       pvtk.PointData(pvtk.Scalars(mat_val_all, name='kernel_value')),
                       'Inversion Grid')
    print '\n\n=================='
    print "WARNING: INDEXING!"
    print '=================='
    vtk.tofile(os.path.join(os.path.curdir, '%s.vtk' % input_file_name.split('_')[0]))

####################### vtk_generator_all #############################


def vtk_generator_all(direname, vertex_file, facet_file):
    """
    VTK file generator out of all the results of all runs
    ATTENTION: ascii.matrix.* should be moved to one dir (direname)
    INFO:
    - facet file is indexed from 0
    - BUT! ascii.* files are indexed from 1, so we need to subtract them (-1)
    - When describing the cells in terms of point indices, the points must be indexed starting at 0.
    """
    print '\n======>> Creating VTK file'
    ascii_files = glob.glob(os.path.join(direname, 'ascii.matrixT.*'))

    print '------> Load Vertex file'
    mesh_points = np.loadtxt(os.path.join(direname, vertex_file), skiprows=2, comments='#')
    print '------> Load Facet file'
    mesh_facets = np.loadtxt(os.path.join(direname, facet_file), dtype=np.int, skiprows=1, comments='#')

    mat_val_all = [0.0]*len(mesh_points)
    for nj in range(len(ascii_files)):
        print '\n------> create VTK file %s' % ascii_files[nj]
        fmatrix = open(os.path.join(ascii_files[nj]), 'r')
        fmatrix_r = fmatrix.readlines()
        mat_indx = []
        mat_val = []
        counter = 0
        for j in range(5, len(fmatrix_r)):
            if counter == 3:
                counter = 0
                continue
            if counter == 0:
                # mat_indx: matrix index
                mat_indx_tmp = fmatrix_r[j].split()
                for i in range(len(mat_indx_tmp)):
                    # IF indexing in ASCII file starts from 0,
                    # we do not need -1.
                    # Otherwise we need it!
                    mat_indx.append(int(mat_indx_tmp[i])-1)
            if counter == 1:
                # mat_val: matrix value
                mat_val_tmp = fmatrix_r[j].split()
                for i in range(len(mat_val_tmp)):
                    mat_val.append(abs(float(mat_val_tmp[i])))
            counter += 1

        for i in range(len(mat_indx)):
            mat_val_all[mat_indx[i]] += mat_val[i]

    vtk = pvtk.VtkData(pvtk.UnstructuredGrid(mesh_points, tetra=mesh_facets),
                       pvtk.PointData(pvtk.Scalars(mat_val_all, name='kernel_value')),
                       'Inversion Grid')
    print '\n\n=================='
    print "WARNING: INDEXING!"
    print '=================='
    vtk.tofile(os.path.join(direname, 'global_all.vtk'))










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
