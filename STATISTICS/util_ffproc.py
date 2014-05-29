import matplotlib.pyplot as plt
import numpy as np
import os

#### XXXX IT HAS A READER: BE CAREFUL ABOUT XCORR AND MEDIAN!!!!

#----------------------reader---------------------------------


def reader(evadd, bands, band_period, all_stations=False, just_high_cc=False, remove_GSN_median=False):
    """
    This function reads the ffproc.ampstt.band....
    """

    GSN_stations = \
        ['II.AAK', 'II.ABKT', 'II.ABPO', 'IU.ADK', 'IU.AFI', 'II.ALE', 'IU.ANMO', 'IU.ANTO', 'CU.ANWB', 'II.ARU',
         'II.ASCN', 'CU.BBGH', 'IU.BBSR', 'CU.BCIP', 'GT.BDFB', 'II.BFO', 'GT.BGCA', 'IU.BILL', 'IC.BJT', 'II.BORG',
         'GT.BOSA', 'II.BRVK', 'IU.CASY', 'IU.CCM', 'IU.CHTO', 'II.CMLA', 'II.COCO', 'IU.COLA', 'IU.COR',
         'GT.CPUP', 'IU.CTAO', 'IU.DAV',  'GT.DBIC', 'II.DGAR', 'IU.DWPF', 'II.EFI', 'IC.ENH', 'II.ERM', 'II.ESK',
         'II.FFC', 'IU.FUNA', 'IU.FURI', 'IU.GNI', 'IU.GRFO', 'CU.GRGR', 'CU.GRTK', 'CU.GTBY', 'IU.GUMO', 'IC.HIA',
         'IU.HKT', 'IU.HNR', 'II.HOPE', 'IU.HRV', 'IU.INCN', 'IU.JOHN', 'II.JTS', 'II.KAPI', 'IU.KBL', 'IU.KBS',
         'II.KDAK', 'IU.KEV', 'IU.KIEV', 'IU.KIP', 'II.KIV', 'IU.KMBO', 'IC.KMI', 'IU.KNTN', 'IU.KONO', 'IU.KOWA',
         'II.KURK', 'II.KWAJ', 'GT.LBTB', 'IU.LCO', 'GT.LPAZ', 'IC.LSA', 'IU.LSZ', 'IU.LVC', 'II.LVZ', 'IU.MA2',
         'IU.MACI', 'IU.MAJO', 'IU.MAKZ', 'II.MBAR', 'IU.MBWA', 'IC.MDJ', 'IU.MIDW', 'II.MSEY', 'IU.MSKU',
         'II.MSVF', 'CU.MTDJ', 'II.NIL', 'II.NNA', 'II.NRIL', 'IU.NWAO', 'II.OBN', 'IU.OTAV', 'IU.PAB', 'II.PALK',
         'IU.PAYG', 'IU.PET', 'II.PFO', 'GT.PLCA', 'IU.PMG', 'IU.PMSA', 'IU.POHA', 'IU.PTCN', 'IU.PTGA', 'IC.QIZ',
         'IU.QSPA', 'IU.RAO', 'IU.RAR', 'II.RAYN', 'IU.RCBR', 'II.RPN', 'IU.RSSD', 'II.SACV', 'IU.SAML',  'IU.SBA',
         'CU.SDDR', 'IU.SDV', 'IU.SFJD', 'II.SHEL', 'IU.SJG', 'IU.SLBS', 'IU.SNZO', 'IC.SSE', 'IU.SSPA', 'II.SUR',
         'IU.TARA', 'IU.TATO', 'II.TAU', 'IU.TEIG', 'CU.TGUH', 'IU.TIXI', 'II.TLY', 'IU.TRIS', 'IU.TRQA', 'IU.TSUM',
         'IU.TUC', 'IU.ULN', 'GT.VNDA', 'IU.WAKE', 'IU.WCI', 'IC.WMQ', 'II.WRAB', 'IU.WVT', 'IC.XAN', 'IU.XMAS',
         'IU.YAK', 'IU.YSS']

    failed = 0
    passed_staev = []
    evlats = []
    evlons = []

    for i in bands:
        try:
            str_i = 'band0' + str(i)
            all_dt_high_cc = []
            dt_GSN = []
            all_dt_event = np.array([])
            all_da_event = np.array([])
            passed_staev_tmp = []

            fio_source = open(os.path.join(evadd, 'outfiles', 'ffproc.source'), 'r')
            f_source = fio_source.readlines()
            ev_year, ev_julianday, ev_hr, ev_min, ev_sec, ev_msec = f_source[1].split()
            evlat, evlon, catalog_depth, inverted_depth = f_source[3].split()
            try:
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[13].split()
            except Exception, e:
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[7].split()

            fio_dt = open(os.path.join(evadd, 'outfiles', 'ffproc.ampstt.' + str_i), 'r')
            f_dt = fio_dt.readlines()
            for j in range(0, len(f_dt)):
                if f_dt[j].strip().startswith('#'):
                    continue
                info_dt = f_dt[j].split()
                lat = float(info_dt[2])
                lon = float(info_dt[3])
                epi = float(info_dt[4])
                sta_id = info_dt[5]
                xcorr = float(info_dt[6])
                dt = float(info_dt[8])
                da = float(info_dt[12])
                clip_taumax = int(info_dt[19])
                # First we collect all the information and then filter it in the next step!
                passed_staev_tmp.append([lat, lon, xcorr, band_period[str(i)], epi, sta_id, clip_taumax])
                all_dt_event = np.append(all_dt_event, dt)
                all_da_event = np.append(all_da_event, da/1.e9)
                if just_high_cc:
                    if xcorr >= just_high_cc:
                        if clip_taumax != 1:
                            all_dt_high_cc.append(dt)
                            if remove_GSN_median:
                                station_id = '%s.%s' % (info_dt[5].split('.')[0], info_dt[5].split('.')[1])
                                if station_id in GSN_stations:
                                    dt_GSN.append(dt)
            if remove_GSN_median and len(dt_GSN) > 10:
                np_median = np.median(dt_GSN)
            elif just_high_cc and len(all_dt_high_cc) > 0:
                np_median = np.median(all_dt_high_cc)
            else:
                np_median = np.median(all_dt_event)

            #np_median = 0
            print 'MEDIAN: %s ... length GSN: %s' % (np_median, len(dt_GSN))

            all_dt_median = all_dt_event - np_median
            all_da_median = all_da_event - np.median(all_da_event)

            for k in range(len(all_dt_median)):
                passed_staev_tmp[k].insert(2, all_dt_median[k])
            for k in range(len(all_da_median)):
                passed_staev_tmp[k].insert(3, all_da_median[k])
            passed_staev.append(passed_staev_tmp)
        except Exception, e:
            print 'ERROR: %s' % e
            failed += 1
    return passed_staev

#----------------------filters---------------------------------


def filters(all_staev, bands, xcorr_limit=False, all_stations=True):
    """
    filters all station-event pairs based on inputs
    """

    GSN_stations = \
        ['II.AAK', 'II.ABKT', 'II.ABPO', 'IU.ADK', 'IU.AFI', 'II.ALE', 'IU.ANMO', 'IU.ANTO', 'CU.ANWB', 'II.ARU',
         'II.ASCN', 'CU.BBGH', 'IU.BBSR', 'CU.BCIP', 'GT.BDFB', 'II.BFO', 'GT.BGCA', 'IU.BILL', 'IC.BJT', 'II.BORG',
         'GT.BOSA', 'II.BRVK', 'IU.CASY', 'IU.CCM', 'IU.CHTO', 'II.CMLA', 'II.COCO', 'IU.COLA', 'IU.COR',
         'GT.CPUP', 'IU.CTAO', 'IU.DAV',  'GT.DBIC', 'II.DGAR', 'IU.DWPF', 'II.EFI', 'IC.ENH', 'II.ERM', 'II.ESK',
         'II.FFC', 'IU.FUNA', 'IU.FURI', 'IU.GNI', 'IU.GRFO', 'CU.GRGR', 'CU.GRTK', 'CU.GTBY', 'IU.GUMO', 'IC.HIA',
         'IU.HKT', 'IU.HNR', 'II.HOPE', 'IU.HRV', 'IU.INCN', 'IU.JOHN', 'II.JTS', 'II.KAPI', 'IU.KBL', 'IU.KBS',
         'II.KDAK', 'IU.KEV', 'IU.KIEV', 'IU.KIP', 'II.KIV', 'IU.KMBO', 'IC.KMI', 'IU.KNTN', 'IU.KONO', 'IU.KOWA',
         'II.KURK', 'II.KWAJ', 'GT.LBTB', 'IU.LCO', 'GT.LPAZ', 'IC.LSA', 'IU.LSZ', 'IU.LVC', 'II.LVZ', 'IU.MA2',
         'IU.MACI', 'IU.MAJO', 'IU.MAKZ', 'II.MBAR', 'IU.MBWA', 'IC.MDJ', 'IU.MIDW', 'II.MSEY', 'IU.MSKU',
         'II.MSVF', 'CU.MTDJ', 'II.NIL', 'II.NNA', 'II.NRIL', 'IU.NWAO', 'II.OBN', 'IU.OTAV', 'IU.PAB', 'II.PALK',
         'IU.PAYG', 'IU.PET', 'II.PFO', 'GT.PLCA', 'IU.PMG', 'IU.PMSA', 'IU.POHA', 'IU.PTCN', 'IU.PTGA', 'IC.QIZ',
         'IU.QSPA', 'IU.RAO', 'IU.RAR', 'II.RAYN', 'IU.RCBR', 'II.RPN', 'IU.RSSD', 'II.SACV', 'IU.SAML',  'IU.SBA',
         'CU.SDDR', 'IU.SDV', 'IU.SFJD', 'II.SHEL', 'IU.SJG', 'IU.SLBS', 'IU.SNZO', 'IC.SSE', 'IU.SSPA', 'II.SUR',
         'IU.TARA', 'IU.TATO', 'II.TAU', 'IU.TEIG', 'CU.TGUH', 'IU.TIXI', 'II.TLY', 'IU.TRIS', 'IU.TRQA', 'IU.TSUM',
         'IU.TUC', 'IU.ULN', 'GT.VNDA', 'IU.WAKE', 'IU.WCI', 'IC.WMQ', 'II.WRAB', 'IU.WVT', 'IC.XAN', 'IU.XMAS',
         'IU.YAK', 'IU.YSS']

    indx_failed = []
    for i in range(len(all_staev[0])):
        flag = True
        # required stations are those that can pass the xcorr test in all bands
        # therefore, it first checks the xcorr in all bands for one station
        # Obviously, bands can be in the shape of [band] in order to check just one band and not all!
        for j in range(len(bands)):
            if all_staev[j][i][4] < xcorr_limit:
                flag = False
                break
            if all_staev[j][i][8] == 1:
                flag = False
                break
            #station_id = '%s.%s' % (all_staev[j][i][7].split('.')[0], all_staev[j][i][7].split('.')[1])
            #if not station_id in GSN_stations:
            #    print '%s not in GSN stations' % station_id
            #    flag = False
            #    break
        if not flag:
            indx_failed.append(i)

    indx_failed.reverse()
    for j in range(len(bands)):
        for i in indx_failed:
            all_staev[j].pop(i)
    return all_staev 

#----------------------ffpscatter---------------------------------


def ffpscatter(passed_staev, all_events=False):
    """
    Plot dT and observed A for all the stations in one event
    """
    if not all_events:
        for i in range(len(passed_staev)):
            for j in range(len(passed_staev[i])):
                plt.plot(passed_staev[i][j][5], passed_staev[i][j][2], 'ko')
        plt.show()
    else:
        ls_all_passed_band = []
        ls_all_passed_dt = []
        for k in range(len(passed_staev)):
            for i in range(len(passed_staev[k])):
                for j in range(len(passed_staev[k][i])):
                    plt.plot(passed_staev[k][i][j][5], passed_staev[k][i][j][2], 'ko')
                    ls_all_passed_band.append(passed_staev[k][i][j][5])
                    ls_all_passed_dt.append(passed_staev[k][i][j][2])

        plt.ion()
        import py2mat_mod
        py2mat_mod.py2mat(ls_all_passed_band, 'dispersion_all_passed_band', 'dispersion_all_passed_band')
        py2mat_mod.py2mat(ls_all_passed_dt, 'dispersion_all_passed_dt', 'dispersion_all_passed_dt')

        plt.show()

#----------------------stamean---------------------------------


def stamean(passed_staev):
    """
    return dt mean for different bands in one event
    """
    per = []
    dt_mean = []
    a_mean = []
    flag = True
    for i in range(len(passed_staev)):
        dt = np.array([])
        a = np.array([])
        for j in range(len(passed_staev[i])):
            dt = np.append(dt, passed_staev[i][j][2])
            a = np.append(a, passed_staev[i][j][3])
        try:
            per.append(passed_staev[i][0][5])
            dt_mean.append(np.mean(dt))
            a_mean.append(np.mean(a))
        except Exception, e:
            print 'ERROR in stamean: %s' % e
            flag = False
    return per, dt_mean, a_mean, flag

#----------------------ffplotmean---------------------------------


def ffplotmean(per, dt_mean):
    """
    Plot mean dT and observed A for all the stations in one event
    """
    plt.plot(per, dt_mean, 'k')
    plt.show()

#----------------------allffplotmean---------------------------------


def allffplotmean(per, all_dt_mean):
    """
    Plot mean dT and observed A for all the stations for each events
    """
    for i in range(len(all_dt_mean)):
        plt.plot(per, all_dt_mean[i][0])
    plt.show()

#----------------------meanall_ffplot---------------------------------


def meanall_ffplot(per, all_dt_mean, all_a_mean, all_tt_single):
    """
    Plot mean dT and observed A for all the stations in all events
    """
    meanall_dt = []
    meanall_a = []
    for j in range(len(all_dt_mean[0][0])):
        weight_tt = 0
        weight_aa = 0
        tt = 0
        aa = 0
        for i in range(len(all_dt_mean)):
            tt += all_dt_mean[i][0][j] * all_dt_mean[i][1]
            weight_tt += all_dt_mean[i][1]
            # some strange observed A exist...
            # this is a very simple way and naive to go around it!
            if not all_a_mean[i][0][j] > 10:
                aa += all_a_mean[i][0][j] * all_a_mean[i][1]
                weight_aa += all_a_mean[i][1]
            else:
                print 'Error: %s for observed A' % all_a_mean[i][0][j]
        meanall_dt.append(tt/weight_tt)
        #meanall_a.append(aa/weight_aa)
    
    ### A TEST
    plt.ion()
    plt.figure()
    plt.subplot(1, 1, 1)
    all_tt_mean = []
    all_tt_std = []
    for i in range(len(all_tt_single)):
        all_tt_mean.append(np.mean(all_tt_single[i]))
        all_tt_std.append(np.std(all_tt_single[i]))
     
    plt.xlabel('Dominant Period', fontsize='x-large', weight='bold')
    x = [2.7, 3.7, 5.3, 7.5, 10.6, 15.0, 21.2, 30.0]
    plt.xlim(xmin=0.0)
    plt.vlines(x, 0.0, 1.6, linestyle='--')
    # !!! change these according to your case!
    plt.xlim(1.7, 31)
    plt.ylim(0.0, 1.6)
    plt.xticks(x, fontsize='x-large', weight='bold')
    plt.yticks(fontsize='x-large', weight='bold')
    plt.plot(per, all_tt_mean, lw=3, color='black', label='Mean')
    plt.plot(per, all_tt_std, lw=3, color='red', label='STD')
    plt.legend(loc=5, prop={'size': 32})
    plt.show()
    ### FINISH A TEST
    
    plt.figure() 
    plt.subplot(1, 1, 1)
    plt.ylabel('Time difference (dT)', fontsize='x-large', weight='bold')
    plt.xlabel('Dominant Period', fontsize='x-large', weight='bold')
    x = [2.7, 3.7, 5.3, 7.5, 10.6, 15.0, 21.2, 30.0]
    plt.xlim(xmin=0.0)
    #plt.ylim(ymin=0.25, ymax=0.65)
    plt.vlines(x, 0.0, 0.65, linestyle='--')
    plt.ylim(0.0, 0.65)
    plt.xlim(1.7, 31.)
    plt.xticks(x, fontsize='x-large', weight='bold')
    plt.yticks(fontsize='x-large', weight='bold')
    #pltitle = evname
    #pltitle += '\nxcorr >= %s' %(xcorr_limit)
    #pltitle = '#station-event pairs (dT): %s\n' %(weight_tt)
    #pltitle += '#station-event pairs (A): %s' %(weight_aa)
    pltitle = '#station-event pairs (dT): %s' % weight_tt
    plt.title(pltitle, fontsize='x-large', weight='bold')
    
    # writing meanall_dt for later usages!
    import pickle
    meanall_file = open(os.path.join('.', 'statistics', 'meanall_dt'), 'w')
    print "start pickling the meanall_file in %s..." % meanall_file
    pickle.dump(meanall_dt, meanall_file)
    print "DONE"

    print '\n\n=========='
    print 'bands:'
    print per
    print 'mean values:'
    print meanall_dt

    import py2mat_mod
    py2mat_mod.py2mat(per, 'dispersion_period', 'dispersion_period')
    py2mat_mod.py2mat(meanall_dt, 'dispersion_meanall_dt', 'dispersion_meanall_dt')

    plt.plot(per, meanall_dt, linewidth=3)

    #plt.subplot(2, 1, 2)
    #plt.ylabel('Observed Amplitude', fontsize = 'x-large', weight = 'bold')
    #plt.xlabel('Dominant Period', fontsize = 'x-large', weight = 'bold')
    #plt.xticks(fontsize = 'x-large', weight = 'bold')
    #plt.yticks(fontsize = 'x-large', weight = 'bold')
    #plt.plot(per, meanall_a, linewidth=3)
    
    plt.show()


#----------------------writer---------------------------------


def writer(passed_staev, bands):
    """
    write the lat, lon, dT, A and xcorr for each station in  different bands
    """
    for i in range(len(passed_staev)):
        try:
            fio = open(os.path.join('.', 'statistics', 'bands0' + str(bands[i])), 'a+')
            for j in range(len(passed_staev[i])):
                msg = '%s,%s,%s,%s,%s,%s\n' %(passed_staev[i][j][0], passed_staev[i][j][1], passed_staev[i][j][2],
                                              passed_staev[i][j][3], passed_staev[i][j][4], len(passed_staev[i]))
                fio.writelines(msg)
            fio.close()
        except Exception, e:
            print 'ERROR in writer %s' % e

# ------------------- round_to --------------------------


def round_to(n, precision):
    """
    rounding the numbers!
    """
    correction = 0.5 if n >= 0 else -0.5
    rounded = int(n/precision+correction)*precision
    rounded2 = round(rounded, 6)
    return rounded2