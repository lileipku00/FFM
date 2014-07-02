#!/usr/bin/env python
# -*- coding: utf-8 -*-

#### XXXX IT HAS A READER: BE CAREFUL ABOUT XCORR AND MEDIAN!!!!

#-------------------------------------------------------------------
#   Filename:  Pdiff_dt_density.py
#   Purpose:   plot density map of all measured dt by FF
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import glob
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import multiprocessing
import numpy as np
from obspy.core.util import locations2degrees
import os
import pickle
import sys
import time

# ------------------- INPUT -----------------------------
#processed_events_add = '/import/neptun-dunkles/hosseini/PROCESSING/ECORR'
#processed_events_add = '/home/hosseini/Work/Scripts/gitHUB/MEASUREMENTS/Pdiff_measure_1_sec_LAMBDA_1-5_90_180'
processed_events_add = '/import/neptun-radler/hosseini-downloads/KASRA/SCRIPTS/gitHUB/myrepo_gitHUB/CORR_MEASUREMENTS'
band = sys.argv[1]
#band = 'BB'
xcorr_limit = 0.8
gr_x = 720
npts = 1800
parts = 30
#gr_x = 18
#npts = 180
#parts = 1
projection = 'robin'
ray_coverage = sys.argv[2]
if ray_coverage == 'F':
    ray_coverage = False
elif ray_coverage == 'T':
    ray_coverage = True
print 'Ray coverage: %s' % ray_coverage
read_only = False

remove_GSN_median = True

# MAP projection
long_0 = 0.

GSN_stations = \
    ['II.AAK', 'II.ABKT', 'II.ABPO', 'IU.ADK', 'IU.AFI', 'II.ALE', 'IU.ANMO', 'IU.ANTO', 'CU.ANWB', 'II.ARU', 'II.ASCN',
        'CU.BBGH', 'IU.BBSR', 'CU.BCIP', 'GT.BDFB', 'II.BFO', 'GT.BGCA', 'IU.BILL', 'IC.BJT', 'II.BORG', 'GT.BOSA',
        'II.BRVK',  'IU.CASY', 'IU.CCM', 'IU.CHTO', 'II.CMLA', 'II.COCO', 'IU.COLA', 'IU.COR', 'GT.CPUP', 'IU.CTAO',
        'IU.DAV',  'GT.DBIC', 'II.DGAR', 'IU.DWPF', 'II.EFI', 'IC.ENH', 'II.ERM', 'II.ESK', 'II.FFC', 'IU.FUNA',
        'IU.FURI', 'IU.GNI', 'IU.GRFO', 'CU.GRGR', 'CU.GRTK', 'CU.GTBY', 'IU.GUMO', 'IC.HIA', 'IU.HKT', 'IU.HNR',
        'II.HOPE', 'IU.HRV', 'IU.INCN', 'IU.JOHN', 'II.JTS', 'II.KAPI', 'IU.KBL', 'IU.KBS', 'II.KDAK', 'IU.KEV',
        'IU.KIEV', 'IU.KIP', 'II.KIV', 'IU.KMBO', 'IC.KMI', 'IU.KNTN', 'IU.KONO', 'IU.KOWA', 'II.KURK', 'II.KWAJ',
        'GT.LBTB', 'IU.LCO', 'GT.LPAZ', 'IC.LSA', 'IU.LSZ', 'IU.LVC', 'II.LVZ', 'IU.MA2', 'IU.MACI', 'IU.MAJO',
        'IU.MAKZ', 'II.MBAR', 'IU.MBWA', 'IC.MDJ', 'IU.MIDW', 'II.MSEY', 'IU.MSKU', 'II.MSVF', 'CU.MTDJ', 'II.NIL',
        'II.NNA', 'II.NRIL', 'IU.NWAO', 'II.OBN', 'IU.OTAV', 'IU.PAB', 'II.PALK', 'IU.PAYG', 'IU.PET', 'II.PFO',
        'GT.PLCA', 'IU.PMG', 'IU.PMSA', 'IU.POHA', 'IU.PTCN', 'IU.PTGA', 'IC.QIZ', 'IU.QSPA', 'IU.RAO', 'IU.RAR',
        'II.RAYN', 'IU.RCBR', 'II.RPN', 'IU.RSSD', 'II.SACV', 'IU.SAML',  'IU.SBA', 'CU.SDDR', 'IU.SDV', 'IU.SFJD',
        'II.SHEL', 'IU.SJG', 'IU.SLBS', 'IU.SNZO', 'IC.SSE', 'IU.SSPA', 'II.SUR', 'IU.TARA', 'IU.TATO', 'II.TAU',
        'IU.TEIG', 'CU.TGUH', 'IU.TIXI', 'II.TLY', 'IU.TRIS', 'IU.TRQA', 'IU.TSUM', 'IU.TUC', 'IU.ULN', 'GT.VNDA',
        'IU.WAKE', 'IU.WCI', 'IC.WMQ', 'II.WRAB', 'IU.WVT', 'IC.XAN', 'IU.XMAS', 'IU.YAK', 'IU.YSS']


# -------------------------------------------------------

#if raw_input('Moved the items in MAP_* dir?(y/n)').lower() == 'n':
#    sys.exit()
#if not read_only:
#    if raw_input('Removed the items in MAP_OUTPUT dir?(y/n)').lower() == 'n':
#        sys.exit()

#############################################################
####################### FUNCTIONS ###########################
#############################################################

# -----------_get_colormap--------------------


def _get_colormap(colors, colormap_name):
    """
    A simple helper function facilitating linear colormap creation.
    """
    # Sort and normalize from 0 to 1.
    indices = np.array(sorted(colors.iterkeys()))
    normalized_indices = (indices - indices.min()) / indices.ptp()

    # Create the colormap dictionary and return the colormap.
    cmap_dict = {"red": [], "green": [], "blue": []}
    for _i, index in enumerate(indices):
        color = colors[index]
        cmap_dict["red"].append((normalized_indices[_i], color[0], color[0]))
        cmap_dict["green"].append((normalized_indices[_i], color[1], color[1]))
        cmap_dict["blue"].append((normalized_indices[_i], color[2], color[2]))
    return LinearSegmentedColormap(colormap_name, cmap_dict)

# -----------point_finder--------------------


def point_finder(lon, lat, grd):
    """
    find the index for a given pair of lon,lat
    --------------------------
    -                        -
    ------> i                -
    -       *                -
    -       j                -
    -       ^                -
    -       |                -
    --------------------------
    """

    for i in range(len(grd[2][0])-1):
        if grd[2][0][i] <= lon < grd[2][0][i+1]:
            break
    for j in range(len(grd[3])-1):
        if grd[3][j][0] <= lat < grd[3][j+1][0]:
            break
    return i, j

# -----------ray_density--------------------


def ray_density(lat1, lon1, lat2, lon2, dt=1, gr_x=360, gr_y=180, npts=180, projection='robin', ray_coverage=False):
    """
    Create the DATA array which contains the info for ray density

    Procedure:
    1. make a grid based on the inputs (grd)
    grd: lon, lat, x, y
    2. find the great circle points:
    note that lon , lat are actually x and y!
    3. calculate the distance and find the middle point
    4. subtracting 97 degrees from the distance and find all the points on that section
    5. data ---> zero array with x*y elements

    """

    global long_0

    mymap = Basemap(projection=projection, lon_0=long_0, lat_0=0)
    #npts=max(gr_x, gr_y)
    # grd[2]: longitude
    # grd[3]: latitude
    grd = mymap.makegrid(gr_x, gr_y, returnxy=True)

    lons, lats = mymap.gcpoints(lon1, lat1, lon2, lat2, npts)
    dist = locations2degrees(lat1, lon1, lat2, lon2)

    # npts points on dist...how many on (dist-97)!: (dist-97)*npts/dist....but we also need to make it half!
    bap = int((dist - 97.0)*npts/dist)/2
    midlon = len(lons)/2
    midlat = len(lats)/2
    lons = lons[midlon-bap:midlon+1+bap]
    lats = lats[midlat-bap:midlat+1+bap]

    data = np.zeros([len(grd[2]), len(grd[3])])
    if not len(lons) == len(lats):
        sys.exit('ERROR: Lengths longitudes and latitudes are not the same! %s and %s' % (len(lons), len(lats)))

    for i in range(len(lons)):
        xi, yi = point_finder(lons[i], lats[i], grd)
        # first one is latitude and second longitude
        try:
            #data[yi][xi] = dt/float(dist-97.0)
            data[yi][xi] += dt/len(lons)
        except Exception, e:
            print '\nException: %s' % e

    if ray_coverage:
        data[np.nonzero(data)] = 1
    return data

# -----------calculator--------------------


def calculator(DATA, passed_staev, gr_x, npts, start, end, projection='robin', ray_coverage=False):
    """
    Calculate the matrix that will be plotted later!
    Matrix can be for ray coverage or for travel-time map!
    """
    nonzero = []
    # Assume the same number if elements in both x and y
    gr_y = gr_x

    for i in range(start, end):
        print i,
        sys.stdout.flush()
        data = ray_density(passed_staev[i][4], passed_staev[i][5], passed_staev[i][0], passed_staev[i][1],
                           dt=passed_staev[i][2], gr_x=gr_x, gr_y=gr_y, npts=npts, projection=projection,
                           ray_coverage=ray_coverage)
        if DATA is None:
            DATA = data.copy()
        else:
            DATA += data
        nonzero_tmp = np.nonzero(data)
        for j in range(len(nonzero_tmp[0])):
            nonzero.append((nonzero_tmp[0][j], nonzero_tmp[1][j]))
    fi = open('MAP_OUTPUT/DATA-%s' % start, 'w')
    pickle.dump(DATA, fi)
    fi.close()

    fi = open('MAP_OUTPUT/nonzero-%s' % start, 'w')
    pickle.dump(nonzero, fi)
    fi.close()

#############################################################
###################### MAIN PROGRAM #########################
#############################################################
if not read_only:
    band_period = {
        'band01': '30.0sec',
        'band02': '21.2sec',
        'band03': '15.0sec',
        'band04': '10.6sec',
        'band05': '7.5sec',
        'band06': '5.3sec',
        'band07': '3.7sec',
        'band08': '2.7sec'}

    if band in band_period:
        band_p = band_period[band]
    elif band == 'BB':
        band_p = 'BB'
    else:
        sys.exit('Wrong band is specified! (check the INPUT)')

    proc_ev_ls = glob.glob(os.path.join(processed_events_add, '*.*.*.*'))
    print '%s processed events found!' % len(proc_ev_ls)

    failed = 0
    passed_staev = []
    evlats = []
    evlons = []

    for i in range(len(proc_ev_ls)):
        evnt = proc_ev_ls[i]
        all_dt_high_cc = []
        dt_GSN = []
        all_dt_event = []
        passed_staev_tmp = []
        try:
            fio_source = open(os.path.join(evnt, 'outfiles', 'ffproc.source'), 'r')
            f_source = fio_source.readlines()
            ev_year, ev_julianday, ev_hr, ev_min, ev_sec, ev_msec = f_source[1].split()
            evlat, evlon, catalog_depth, inverted_depth = f_source[3].split()
            try:
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[13].split()
            except Exception, e:
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[7].split()

            fio_dt = open(os.path.join(evnt, 'outfiles', 'ffproc.ampstt.' + band), 'r')
            f_dt = fio_dt.readlines()
            for j in range(0, len(f_dt)):
                if f_dt[j].strip().startswith('#'):
                    continue
                info_dt = f_dt[j].split()
                lat = float(info_dt[2])
                lon = float(info_dt[3])
                xcorr = float(info_dt[6])
                dt = float(info_dt[8])
                clip_taumax = int(info_dt[19])
                all_dt_event.append(dt)
                if xcorr >= xcorr_limit:
                    if clip_taumax != 1:
                        passed_staev_tmp.append([lat, lon, xcorr, float(evlat), float(evlon), i])
                        all_dt_high_cc.append(dt)
                        if remove_GSN_median:
                            station_id = '%s.%s' % (info_dt[5].split('.')[0], info_dt[5].split('.')[1])
                            if station_id in GSN_stations:
                                dt_GSN.append(dt)
            if remove_GSN_median and len(dt_GSN) > 10:
                np_median = np.median(dt_GSN)
            elif len(all_dt_high_cc) > 0:
                np_median = np.median(all_dt_high_cc)
            else:
                np_median = np.median(all_dt_event)

            all_dt_median = all_dt_high_cc - np_median
            for k in range(len(all_dt_median)):
                passed_staev_tmp[k].insert(2, all_dt_median[k])
                passed_staev.append(passed_staev_tmp[k])
            evlats.append(float(evlat))
            evlons.append(float(evlon))
        except Exception, e:
            print 'ERROR: %s' % e
            failed += 1

    xcorr_sum = 0
    dt_sum = 0
    for i in range(len(passed_staev)):
        xcorr_sum += passed_staev[i][3]
        dt_sum += passed_staev[i][2]

    print '\n------------------------'
    print '%s events failed!' % failed
    print '%s station-event pairs found...' % len(passed_staev)
    print 'mean xcorr: %s' % (xcorr_sum/len(passed_staev))
    print 'mean dt: %s' % (dt_sum/len(passed_staev))
    print '------------------------'

    #passed_staev = passed_staev[0:1]
    start = 0
    end = len(passed_staev)
    step = (end - start) / parts + 1

    jobs = []
    DATA = None
    for index in xrange(parts):
        starti = start+index*step
        endi = min(start+(index+1)*step, end)
        p = multiprocessing.Process(target=calculator, args=(DATA, passed_staev, gr_x, npts, starti, endi,
                                                             projection, ray_coverage))
        jobs.append(p)
    for i in range(len(jobs)):
        jobs[i].start()

    pp_flag = True
    while pp_flag:
        for proc in jobs:
            if proc.is_alive():
                print '\nRunning...proc: %s' % proc
                time.sleep(10)
                pp_flag = True
                break
            else:
                pp_flag = False
        if not pp_flag:
            print '\nAll the processes are finished...'

# =================== colormap ==============================
# A pretty colormap for use in tomography.
tomo_colormap = _get_colormap({
    0.0: [0.1, 0.0, 0.0], # Reddish black
    0.2: [0.8, 0.0, 0.0],
    0.3: [1.0, 0.7, 0.0],
    0.48: [0.92, 0.92, 0.92],
    0.5: [0.92, 0.92, 0.92], # Light gray
    0.52: [0.92, 0.92, 0.92],
    0.7: [0.0, 0.6, 0.7],
    0.8: [0.0, 0.0, 0.8],
    1.0: [0.0, 0.0, 0.1]},
    "seismic_tomography") # Blueish black
tomo_colormap_2 = _get_colormap({
        0.00000: [0.12941, 0.40000, 0.67451],
        0.11110: [0.12941, 0.40000, 0.67451],
        0.11110: [0.26275, 0.57647, 0.76471],
        0.22220: [0.26275, 0.57647, 0.76471],
        0.22220: [0.57255, 0.77255, 0.87059],
        0.33330: [0.57255, 0.77255, 0.87059],
        0.33330: [0.81961, 0.89804, 0.94118],
        0.44440: [0.81961, 0.89804, 0.94118],
        0.44440: [0.96863, 0.96863, 0.96863],
        0.55560: [1.0, 1.0, 1.0],
        0.55560: [1.0, 1.0, 1.0],
        0.66670: [0.99216, 0.85882, 0.78039],
        0.66670: [0.95686, 0.64706, 0.50980],
        0.77780: [0.95686, 0.64706, 0.50980],
        0.77780: [0.83922, 0.37647, 0.30196],
        0.88890: [0.83922, 0.37647, 0.30196],
        0.88890: [0.69804, 0.09412, 0.16863],
        1.00000: [0.69804, 0.09412, 0.16863]},
        "seismic_tomography_2")
# =================== END colormap ==============================

mymap = Basemap(projection=projection, lon_0=long_0, lat_0=0)
#mymap.drawmapboundary(fill_color = 'black', color = 'red')
#mymap.drawcoastlines(color='black')

data_ls = glob.glob('MAP_OUTPUT/DATA-*')
nonzero_ls = glob.glob('MAP_OUTPUT/nonzero-*')
DATA = None
nonzero = None

print '\n\n-------------'
print 'Pickle load DATA and NONZERO...'
for i in range(len(data_ls)):
    print i,
    data = pickle.load(open(data_ls[i]))
    nonzero_tmp = pickle.load(open(nonzero_ls[i]))
    if DATA is None:
        DATA = data
    else:
        DATA += data
    if nonzero is None:
        nonzero = nonzero_tmp
    else:
        for noz_item in nonzero_tmp: 
            nonzero.append(noz_item)

print '-------------'
print 'Averaging over each cell...'
enum = 1
nonzero_unique = []
nonzero.sort()
for i in range(len(nonzero)-1, -1, -1):
    if nonzero[i] == nonzero[i-1]:
        enum += 1
    else:
        nonzero_unique.append([nonzero[i], enum])
        enum = 1

if not ray_coverage:
    for i in range(len(nonzero_unique)):
        DATA[nonzero_unique[i][0]] = DATA[nonzero_unique[i][0]]/float(nonzero_unique[i][1])

print '-------------'
print 'Plotting...'
enum = 1
gr_y = gr_x
grd = mymap.makegrid(gr_x, gr_y, returnxy=True)
#mymap.contourf(grd[2], grd[3], DATA)

### vmin = max(abs(np.min(DATA)), abs(np.max(DATA)))
### if not ray_coverage:
###     mymap.pcolormesh(grd[2], grd[3], DATA, cmap=tomo_colormap_2, vmin=-1*vmin/100., vmax=vmin/100.)
### else:
###     import matplotlib.cm as cm
###     mymap.pcolormesh(grd[2], grd[3], DATA, cmap=cm.gray, vmax=10)

### #plt.hexbin(grd[2], grd[3], DATA)
### cbar = plt.colorbar(orientation='horizontal')
### cbar.ax.tick_params(labelsize=12) 
### plt.show()
### 
### 
### # Clean 1.0% (dT)
### from scipy.ndimage import gaussian_filter
### DATA_filt = gaussian_filter(DATA, sigma=10.0, order=0)
### DATA_filt = DATA
### mymap = Basemap(projection=projection, lon_0=long_0, lat_0=0)
### mymap.drawcoastlines()
### # this one is for CMB, based on your notes and corresponding to 1.%
### mymap.pcolormesh(grd[2], grd[3], DATA_filt, cmap=tomo_colormap_2, vmin=-0.02537847, vmax=0.02537847)
### cbar = plt.colorbar(orientation='horizontal')
### cbar.ax.tick_params(labelsize=24)
### cbar.ax.set_xticklabels(['-1.0%', ' ', ' ', ' ', '0%', ' ', ' ', ' ', '1.0%'])
### plt.show()
### 
### # Clean 1.5%
### from scipy.ndimage import gaussian_filter
### DATA_filt = gaussian_filter(DATA, sigma=5.0, order=0)
### DATA_filt = DATA
### mymap = Basemap(projection=projection, lon_0=long_0, lat_0=0)
### mymap.drawcoastlines()
### # this one is for CMB, based on your notes and corresponding to 1.5%
### mymap.pcolormesh(grd[2], grd[3], DATA_filt, cmap=tomo_colormap_2, vmin=-0.0380677, vmax=0.0380677)
### cbar = plt.colorbar(orientation='horizontal')
### cbar.ax.tick_params(labelsize=16)
### cbar.ax.set_xticklabels(['-1.5%', ' ', ' ', ' ', '0%', ' ', ' ', ' ', '1.5%'])
### plt.show()
### 
### 
### # FOR RAY COVERAGE
### from scipy.ndimage import gaussian_filter
### from matplotlib.colors import LogNorm
### DATA_filt = gaussian_filter(DATA, sigma=5.0, order=0)
### DATA_filt = DATA
### mymap = Basemap(projection=projection, lon_0=long_0, lat_0=0)
### mymap.drawcoastlines(color='black')
### vmin = max(abs(np.min(DATA)), abs(np.max(DATA)))
### print vmin
### mymap.pcolormesh(grd[2], grd[3], DATA_filt, norm=LogNorm(vmin=1, vmax=vmin))
### cbar = plt.colorbar(orientation='horizontal')
### cbar.ax.tick_params(labelsize=24)
### #cbar.ax.set_xticklabels(['0.1', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '285'])
### plt.show()

# ==================== PICKLING THE REQUIRED INFORMATION ===============================
print '-------------'
print "Pickling the final result in the following address:\n%s" % os.path.basename(processed_events_add)
os.mkdir(os.path.join('.', 'MAP_' + os.path.basename(processed_events_add)))
data_fi = open(os.path.join('.', 'MAP_' + os.path.basename(processed_events_add), 'DATA'), 'w')
projec_fi = open(os.path.join('.', 'MAP_' + os.path.basename(processed_events_add), 'projection'), 'w')
grd_fi = open(os.path.join('.', 'MAP_' + os.path.basename(processed_events_add), 'grd'), 'w')
tomo_color_fi = open(os.path.join('.', 'MAP_' + os.path.basename(processed_events_add), 'tomo_colormap_2'), 'w')
long_0_fi = open(os.path.join('.', 'MAP_' + os.path.basename(processed_events_add), 'long_0'), 'w')

pickle.dump(DATA, data_fi)
pickle.dump(projection, projec_fi)
pickle.dump(grd, grd_fi)
pickle.dump(tomo_colormap_2, tomo_color_fi)
pickle.dump(long_0, long_0_fi)

data_fi.close()
projec_fi.close()
grd_fi.close()
tomo_color_fi.close()
long_0_fi.close()

print "DONT FORGET TO MOVE THE %s TO ANOTHER DIRECTORY..." % (os.path.basename(processed_events_add) + '_MAP')

print '''
# For opening the files!
# 24 lines
import glob
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import multiprocessing
import numpy as np
from obspy.core.util import locations2degrees
import os
import pickle
import sys
import time

data_fi = open(os.path.join('.', 'DATA'), 'r')
projec_fi = open(os.path.join('.', 'projection'), 'r')
grd_fi = open(os.path.join('.', 'grd'), 'r')
tomo_color_fi = open(os.path.join('.', 'tomo_colormap_2'), 'r')
long_0_fi = open(os.path.join('.', 'long_0'), 'r')

DATA = pickle.load(data_fi)
projection = pickle.load(projec_fi)
grd = pickle.load(grd_fi)
tomo_colormap_2 = pickle.load(tomo_color_fi)
long_0 = pickle.load(long_0_fi)
'''
# ==================== END PICKLING THE REQUIRED INFORMATION ===============================










# =========================================================================================
# ============================== OTHERS (use with care) ===================================
# =========================================================================================
# import scipy.ndimage as ndimage
# DATA_filt = ndimage.gaussian_filter(DATA, sigma=10.0, order=0)
# mymap = Basemap(projection=projection, lon_0=180, lat_0=0)
# mymap.drawcoastlines()
# mymap.pcolormesh(grd[2], grd[3], -1*DATA_filt, cmap=tomo_colormap, vmin=-1*vmin/100., vmax=vmin/100.)
# plt.colorbar()
# #plt.colorbar(orientation="horizontal")
# plt.show()
# 
# import scipy.ndimage as ndimage
# DATA_filt = ndimage.gaussian_filter(DATA, sigma=5.0, order=0)
# mymap = Basemap(projection=projection, lon_0=180, lat_0=0)
# mymap.drawcoastlines()
# mymap.pcolormesh(grd[2], grd[3], -1*DATA_filt, cmap=tomo_colormap, vmin=-1*vmin/100., vmax=vmin/100.)
# plt.colorbar()
# #plt.colorbar(orientation="horizontal")
# plt.show()

# import scipy.ndimage as ndimage
# DATA_filt = ndimage.gaussian_filter(DATA, sigma=5.0, order=0)
# mymap = Basemap(projection=projection, lon_0=180, lat_0=0)
# mymap.drawcoastlines()
# mymap.pcolormesh(grd[2], grd[3], DATA_filt, cmap=tomo_colormap_2, vmin=-1*vmin/100., vmax=vmin/100.)
# plt.colorbar()
# #plt.colorbar(orientation="horizontal")
# plt.show()
# 
#import scipy.ndimage as ndimage
#from matplotlib.colors import LogNorm
#DATA_filt = ndimage.gaussian_filter(DATA, sigma=5.0, order=0)
#DATA_filt = DATA
#mymap = Basemap(projection=projection, lon_0=180, lat_0=0)
#mymap.drawcoastlines(color='white')

#mymap.pcolormesh(grd[2], grd[3], DATA_filt, cmap=cm.gray, vmin=0, vmax=1.0)
#mymap.pcolormesh(grd[2], grd[3], DATA_filt, norm=LogNorm(vmin=0.1, vmax=285))
#plt.colorbar()
#plt.colorbar(orientation="horizontal")
#cbar.ax.tick_params(labelsize=12) 
#plt.show()

########start = 0
########end = len(passed_staev)
########
######### Divide the task into 64 subtasks
########step = (end - start) / parts + 1
########
######### Create jobserver
########job_server = pp.Server()
########
######### Execute the same task with different amount of active workers and measure the time
########job_server.set_ncpus(ncpus)
########jobs = []
########start_time = time.time()
########print "Starting ", job_server.get_ncpus(), " workers"
########for index in xrange(parts):
########    starti = start+index*step
########    endi = min(start+(index+1)*step, end)
########    # Submit a job which will calculate partial sum 
########    # part_sum - the function
########    # (starti, endi) - tuple with arguments for part_sum
########    # () - tuple with functions on which function part_sum depends
########    # () - tuple with module names which must be imported before part_sum execution
########    print index
########    jobs.append(job_server.submit(calculator, (passed_staev, gr_x, npts, starti, endi), depfuncs=(ray_density,point_finder), modules=('from mpl_toolkits.basemap import Basemap','from obspy.core.util import locations2degrees', 'import numpy as np'), callback=True))
########
#########for job in jobs:
#########   job()
#########job_server.wait()
########print 'I AM HERE!'
######### Retrieve all the results and calculate their sum
########pp_flag = True
########while pp_flag:
########    for proc in jobs:
########        if not proc.finished:
########            print 'At least one is not finished yet'
########            print proc
########            time.sleep(10)
########            pp_flag = True
########            break
########        else:
########            pp_flag = False

### tomo_colormap_2 = _get_colormap({
###     0.00000: [0.12941, 0.40000, 0.67451],  
###     0.11110: [0.12941, 0.40000, 0.67451],
###     0.11110: [0.26275, 0.57647, 0.76471],
###     0.22220: [0.26275, 0.57647, 0.76471],
###     0.22220: [0.57255, 0.77255, 0.87059],
###     0.33330: [0.57255, 0.77255, 0.87059],
###     0.33330: [0.81961, 0.89804, 0.94118],
###     0.44440: [0.81961, 0.89804, 0.94118],
###     0.44440: [0.96863, 0.96863, 0.96863],
###     0.55560: [0.96863, 0.96863, 0.96863],
###     0.55560: [0.99216, 0.85882, 0.78039],
###     0.66670: [0.99216, 0.85882, 0.78039],
###     0.66670: [0.95686, 0.64706, 0.50980],
###     0.77780: [0.95686, 0.64706, 0.50980],
###     0.77780: [0.83922, 0.37647, 0.30196],
###     0.88890: [0.83922, 0.37647, 0.30196],
###     0.88890: [0.69804, 0.09412, 0.16863],
###     1.00000: [0.69804, 0.09412, 0.16863]},
###     "seismic_tomography_2")

#tomo_colormap_2 = _get_colormap({
#    0.00000: [0.69804, 0.09412, 0.16863],  
#    0.11110: [0.69804, 0.09412, 0.16863],
#    0.11110: [0.83922, 0.37647, 0.30196],
#    0.22220: [0.83922, 0.37647, 0.30196],
#    0.22220: [0.95686, 0.64706, 0.50980],
#    0.33330: [0.95686, 0.64706, 0.50980],
#    0.33330: [0.99216, 0.85882, 0.78039],
#    0.44440: [0.99216, 0.85882, 0.78039],
#    0.44440: [0.96863, 0.96863, 0.96863],
#    0.55560: [0.96863, 0.96863, 0.96863],
#    0.55560: [0.81961, 0.89804, 0.94118],
#    0.66670: [0.81961, 0.89804, 0.94118],
#    0.66670: [0.57255, 0.77255, 0.87059],
#    0.77780: [0.57255, 0.77255, 0.87059],
#    0.77780: [0.26275, 0.57647, 0.76471],
#    0.88890: [0.26275, 0.57647, 0.76471],
#    0.88890: [0.12941, 0.40000, 0.67451],
#    1.00000: [0.12941, 0.40000, 0.67451]},
#    "seismic_tomography_2")
