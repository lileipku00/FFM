import matplotlib.pyplot as plt
import numpy as np
import os

#----------------------reader---------------------------------
def reader(evadd, bands, band_period, all_stations=True):
    '''
    This function reads the ffproc.ampstt.band....
    '''
    target_network = ['II', 'IU', 'CU', 'GT', 'IC']
    passed_staev = []
    for i in bands:
        try:
            str_i = 'band0' + str(i)
            passed_staev_tmp = []
            all_dt_event = np.array([])
            all_da_event = np.array([])
            fio_dt = open(os.path.join(evadd, 'outfiles', 
                            'ffproc.ampstt.' + str_i), 'r')
            f_dt = fio_dt.readlines()
            for j in range(2, len(f_dt)):
                info_dt = f_dt[j].split()
                if not all_stations:
                    if not info_dt[9].split('.')[0] in target_network:
                        continue
                xcorr = float(info_dt[2])
                da = float(info_dt[3])
                dt = float(info_dt[5])
                lat = float(info_dt[6])
                lon = float(info_dt[7])
                sta_id = info_dt[9]
                passed_staev_tmp.append([lat, lon, xcorr, band_period[str(i)], sta_id])
                all_dt_event = np.append(all_dt_event, dt)
                all_da_event = np.append(all_da_event, da/1.e9)
            all_dt_median = all_dt_event - np.median(all_dt_event)
            all_da_median = all_da_event - np.median(all_da_event)
            for k in range(len(all_dt_median)):
                passed_staev_tmp[k].insert(2, all_dt_median[k])
            for k in range(len(all_da_median)):
                passed_staev_tmp[k].insert(3, all_da_median[k])
            passed_staev.append(passed_staev_tmp)
        except Exception, e:
            print e
    return passed_staev

#----------------------filters---------------------------------
def filters(all_staev, bands, xcorr_limit = False, all_stations=True):
    '''
    filters all station-event pairs based on priori
    '''
    #target_network = ['II', 'IU', 'CU', 'GT', 'IC']
    ind_failed = []
    for i in range(len(all_staev[0])):
        flag = 'T'
        # required stations are those that can pass the xcorr test in all bands
        # therefore, it first checks the xcorr in all bands for one station
        # then continue calculating the RMS fit to the data
        for j in range(len(bands)):
            if all_staev[j][i][4] < xcorr_limit:
                    flag = 'F'
            #if not all_stations:
            #    if not all_staev[j][i][6].split('.')[0] in target_network:
            #            flag = 'F'
        if flag == 'F':
            ind_failed.append(i)
    ind_failed.reverse()
    for j in range(len(bands)):
        for i in ind_failed:
            all_staev[j].pop(i)
    return all_staev 

#----------------------ffpscatter---------------------------------
def ffpscatter(passed_staev, all_events=False):
    '''
    Plot dT and observed A for all the stations in one event
    '''
    if not all_events:
        for i in range(len(passed_staev)):
            for j in range(len(passed_staev[i])):
                plt.plot(passed_staev[i][j][5], passed_staev[i][j][2], 'ko')
        plt.show()
    else:
        for k in range(len(passed_staev)):
            for i in range(len(passed_staev[k])):
                for j in range(len(passed_staev[k][i])):
                    plt.plot(passed_staev[k][i][j][5], passed_staev[k][i][j][2], 'ko')
        plt.show()

#----------------------stamean---------------------------------
def stamean(passed_staev):
    '''
    return dt mean for different bands in one event
    '''
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
            print e
            flag = False
    return per, dt_mean, a_mean, flag

#----------------------ffplotmean---------------------------------
def ffplotmean(per, dt_mean):
    '''
    Plot mean dT and observed A for all the stations in one event
    '''
    plt.plot(per, dt_mean, 'k')
    plt.show()

#----------------------allffplotmean---------------------------------
def allffplotmean(per, all_dt_mean):
    '''
    Plot mean dT and observed A for all the stations for each events
    '''
    for i in range(len(all_dt_mean)):
        plt.plot(per, all_dt_mean[i][0])
    plt.show()

#----------------------meanall_ffplot---------------------------------
def meanall_ffplot(per, all_dt_mean, all_a_mean):
    '''
    Plot mean dT and observed A for all the stations in all events
    '''
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
            if not all_a_mean[i][0][j]>10:
                aa += all_a_mean[i][0][j] * all_a_mean[i][1]
                weight_aa += all_a_mean[i][1]
            else:
                print 'Error: %s for observed A' %(all_a_mean[i][0][j])
        meanall_dt.append(tt/weight_tt)
        #meanall_a.append(aa/weight_aa)
    
    plt.subplot(1, 1, 1)
    plt.ylabel('Time difference (dT)', fontsize = 'x-large', weight = 'bold')
    plt.xlabel('Dominant Period', fontsize = 'x-large', weight = 'bold')
    x = [2.7, 3.7, 5.3, 7.5, 10.6, 15.0, 21.2, 30.0]
    plt.xlim(xmin=0.0)
    #plt.ylim(ymin=0.25, ymax=0.65)
    plt.vlines(x, 0.0, 0.65)
    plt.xticks(x, fontsize = 'x-large', weight = 'bold')
    plt.yticks(fontsize = 'x-large', weight = 'bold')
    #pltitle = evname
    #pltitle += '\nxcorr >= %s' %(xcorr_limit)
    #pltitle = '#station-event pairs (dT): %s\n' %(weight_tt)
    #pltitle += '#station-event pairs (A): %s' %(weight_aa)
    pltitle = '#station-event pairs (dT): %s' %(weight_tt)
    plt.title(pltitle, fontsize = 'x-large', weight = 'bold')
    
    # writing meanall_dt for later usages!
    import pickle
    meanall_file = open(os.path.join('.', 'statistics', 'meanall_dt'), 'w')
    print "start pickling the meanall_file in %s..." %(meanall_file)
    pickle.dump(meanall_dt, meanall_file)
    print "DONE"
    
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
    '''
    write the lat, lon, dT, A and xcorr for each station in 
    different bands
    '''
    for i in range(len(passed_staev)):
        try:
            fio = open(os.path.join('.', 'statistics', 'bands0' + str(bands[i])), 'a+')
            for j in range(len(passed_staev[i])):
                msg = '%s,%s,%s,%s,%s,%s\n' %(passed_staev[i][j][0], passed_staev[i][j][1],
                                               passed_staev[i][j][2], passed_staev[i][j][3],
                                               passed_staev[i][j][4], len(passed_staev[i]))
                fio.writelines(msg)
            fio.close()
        except Exception, e:
            print e
