import glob
import numpy as np
import os

####################### output_file_reader #############################


def output_file_reader(add_output, req_band='band01'):
    """
    This function reads one output file and pass all the values for further analysis in next steps
    """
    empty_array = np.array([])
    if not os.path.isfile(os.path.join(add_output, 'outfiles', 'ffproc.ampstt.%s' % req_band)):
        print "%s is not found!" % os.path.join(add_output, 'outfiles', 'ffproc.ampstt.%s' % req_band)
        return empty_array
    rd_output = np.loadtxt(os.path.join(add_output, 'outfiles', 'ffproc.ampstt.%s' % req_band), dtype='S', comments='#')
    new_col_cr = np.empty([np.shape(rd_output)[0], 2], dtype=object)
    new_col_cr[:, 0] = os.path.basename(add_output)
    new_col_cr[:, 1] = req_band
    output_sta_evname = np.append(rd_output, new_col_cr, 1)
    return output_sta_evname

####################### event_filter #############################


def event_filter(par_add, min_dp=-10, max_dp=1000):
    """
    Filters the events in one par_dir based on the required inputs
    """
    if not os.path.isdir(par_add):
        print "%s is not a valid directory!" % par_add
        return False, False

    event_adds = glob.glob(os.path.join(par_add, '*.*.*.*'))
    passed_event_adds = []
    passed_event_depths = []
    for i in range(len(event_adds)):
        add_flag = False
        try:
            fio_source = open(os.path.join(event_adds[i], 'outfiles', 'ampinv.source'), 'r')
            f_source = fio_source.readlines()
            ev_year, ev_julianday, ev_hr, ev_min, ev_sec, ev_msec = f_source[1].split()
            evlat, evlon, catalog_depth, inverted_depth = f_source[3].split()
            try:
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[13].split()
            except Exception, e:
                print 'WARNING: inverted moment tensor was not found: %s' % event_adds[i]
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[7].split()
        except Exception, e:
            print 'ERROR: %s' % e

        # Check for depth
        if min_dp <= float(inverted_depth) < max_dp:
            add_flag = True
        if add_flag:
            passed_event_adds.append(event_adds[i])
            passed_event_depths.append(inverted_depth)
    return passed_event_adds, passed_event_depths

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
