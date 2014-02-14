#!/usr/bin/env python
# -*- coding: utf-8 -*-

#### XXXX IT HAS A READER: BE CAREFUL ABOUT XCORR AND MEDIAN!!!!

# ATTENTION: This script is the modified version of P and Pdiff ray_density
# The main goal is to compare these two phases togehter!

#-------------------------------------------------------------------
#   Filename:  P_Pdiff_dt_density.py
#   Purpose:   plot dt calculated with FFM in a density map
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
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import multiprocessing
import numpy as np
from obspy.core.util import locations2degrees
import os
import pickle
import sys
import time

# ------------------- INPUT -----------------------------
processed_events_add = '/import/neptun-helles/hosseini/FFM/P_measure_2_sec_LAMBDA_1-5'
band = 'band01'
#band = 'BB'
xcorr_limit = 0.8
gr_x = 720
npts = 1800
parts = 80
#gr_x = 180
#npts = 1800
#parts = 80
projection = 'robin'
ray_coverage = False
read_only = False

# MAP projection
long_0 = 180
# -------------------------------------------------------

if raw_input('Moved the items in MAP_* dir?(y/n)').lower() == 'n':
    sys.exit()
if not read_only:
    if raw_input('Removed the items in MAP_OUTPUT dir?(y/n)').lower() == 'n':
        sys.exit()

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
    '''
    find the index for a given pair of lon,lat
    '''
    for i in range(len(grd[2][0])-1):
        if grd[2][0][i]<=lon<grd[2][0][i+1]:
            break
    for j in range(len(grd[3])-1):
        if grd[3][j][0]<=lat<grd[3][j+1][0]:
            break
    return i, j

# -----------ray_density--------------------
def ray_density(lat1, lon1, lat2, lon2,
                dt=1, gr_x=360, gr_y=180, 
                npts=180, projection='robin', 
                ray_coverage=False):
    '''
    Create the DATA array which contains the
    info for ray density
    '''
    global long_0
    exist_flag = False

    mymap = Basemap(projection=projection, lon_0=long_0, lat_0=0)
    #npts=max(gr_x, gr_y)
    # grd[2]: longitude
    # grd[3]: latitude
    grd = mymap.makegrid(gr_x, gr_y, returnxy=True)
    data = np.zeros([len(grd[2]), len(grd[3])])

    lons, lats = mymap.gcpoints(lon1, lat1, lon2, lat2, npts)
    dist = locations2degrees(lat1, lon1, lat2, lon2)
    if 80.0<=dist<=85.0:
        print dist,
        exist_flag = True
    else:
        print 'None',
        return data, exist_flag
    bap = int((dist - 70.0)*npts/dist)/2
    midlon = len(lons)/2
    midlat = len(lats)/2
    lons = lons[midlon-bap:midlon+1+bap]
    lats = lats[midlat-bap:midlat+1+bap]
    for i in range(len(lons)):
        xi, yi = point_finder(lons[i], lats[i], grd)
        # first one is latitude and second longitude
        try:
            #data[yi][xi] = dt/float(dist-97.0)
            data[yi][xi] += dt/len(lons)
        except Exception, e:
            print e
    if ray_coverage:
        data[np.nonzero(data)] = 1
    return data, exist_flag

# -----------calculator--------------------
def calculator(DATA, passed_staev, gr_x, npts, start, end, 
                    projection='robin', ray_coverage=False): 
    global long_0

    mymap = Basemap(projection=projection, lon_0=long_0, lat_0=0)
    nonzero = []
    gr_y = gr_x
    grd = mymap.makegrid(gr_x, gr_y, returnxy=True)
    for i in range(start, end):
        print i,
        sys.stdout.flush()
        data, exist_flag = ray_density(passed_staev[i][4], passed_staev[i][5],
                        passed_staev[i][0], passed_staev[i][1], 
                        dt=passed_staev[i][2], 
                        gr_x=gr_x, gr_y=gr_y, npts=npts, 
                        projection=projection, 
                        ray_coverage=ray_coverage)
        if not i == end-1:
            if not exist_flag: 
                continue
        if DATA is None: DATA = data.copy()
        else: DATA += data
        nonzero_tmp = np.nonzero(data)
        for j in range(len(nonzero_tmp[0])):
            nonzero.append((nonzero_tmp[0][j], nonzero_tmp[1][j]))
    fi = open('MAP_OUTPUT/DATA-' + str(start), 'w')
    pickle.dump(DATA, fi)
    fi.close()
    fi = open('MAP_OUTPUT/nonzero-' + str(start), 'w')
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
    print '%s processed events found!' %(len(proc_ev_ls))
    failed = 0
    passed_staev = []
    evlats = []
    evlons = []
    print '\nERRORS:'
    for i in range(len(proc_ev_ls)):
        evnt = proc_ev_ls[i]
        all_dt_high_cc = []
        all_dt_event = []
        passed_staev_tmp = []
        try:
            fio_dt = open(os.path.join(evnt, 'outfiles', 'ffproc.ampstt.' + band), 'r')
            fio_source = open(os.path.join(evnt, 'outfiles', 'ampinv.source'), 'r')
            f_source = fio_source.readlines()
            ev_year, ev_julianday, ev_hr, ev_min, ev_sec, ev_msec = f_source[1].split()
            evlat, evlon, catalog_depth, inverted_depth = f_source[3].split()
            try:
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[13].split()
            except Exception, e:
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[7].split()
            f_dt = fio_dt.readlines()
            for j in range(2, len(f_dt)):
                info_dt = f_dt[j].split()
                xcorr = float(info_dt[2])
                dt = float(info_dt[5])
                lat = float(info_dt[6])
                lon = float(info_dt[7])
                all_dt_event.append(dt)
                if xcorr >= xcorr_limit:
                    passed_staev_tmp.append([lat, lon, xcorr, float(evlat), float(evlon), i])
                    all_dt_high_cc.append(dt)
            if len(all_dt_high_cc) > 0:
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
            print e
            failed += 1

    print '\n------------------------'
    xcorr_sum = 0
    dt_sum = 0
    for i in range(len(passed_staev)):
        xcorr_sum += passed_staev[i][3]
        dt_sum += passed_staev[i][2]

    print '%s events failed!' %(failed)
    print '%s station-event pairs found...' %(len(passed_staev))
    print 'mean xcorr: %s' %(xcorr_sum/len(passed_staev))
    print 'mean dt: %s' %(dt_sum/len(passed_staev))
    print '------------------------'

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

    
    mymap = Basemap(projection=projection, lon_0=long_0, lat_0=0)
    #mymap.drawmapboundary(zorder=100)
    mymap.drawcoastlines()
    #passed_staev = passed_staev[0:100]
    start = 0
    end = len(passed_staev)
    step = (end - start) / parts + 1
    jobs = []
    DATA=None
    for index in xrange(parts):
        starti = start+index*step
        endi = min(start+(index+1)*step, end)
        p = multiprocessing.Process(target=calculator, 
                    args=(DATA, passed_staev, gr_x, npts, starti, 
                            endi, projection, ray_coverage))
        jobs.append(p)
    for i in range(len(jobs)):
        jobs[i].start()

    pp_flag = True
    while pp_flag:
        for proc in jobs:
            if proc.is_alive():
                print '\nRunning...'
                print proc
                time.sleep(10)
                pp_flag = True
                break
            else:
                pp_flag = False


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

mymap = Basemap(projection=projection, lon_0=long_0, lat_0=0)
#mymap.drawmapboundary(fill_color = 'black', color = 'red')
mymap.drawcoastlines(color='black')

data_ls = glob.glob('MAP_OUTPUT/DATA-*')
nonzero_ls = glob.glob('MAP_OUTPUT/nonzero-*')
DATA = None
nonzero = None
print '-------------'
print '\n\nPickle...'
for i in range(len(data_ls)):
    print i
    data = pickle.load(open(data_ls[i]))
    nonzero_tmp = pickle.load(open(nonzero_ls[i]))
    if DATA is None: DATA = data
    else: DATA += data 
    if nonzero is None: nonzero = nonzero_tmp
    else:
        for noz_item in nonzero_tmp: 
            nonzero.append(noz_item)

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
        DATA[nonzero_unique[i][0]] = DATA[nonzero_unique[i][0]]/nonzero_unique[i][1]

gr_y = gr_x
grd = mymap.makegrid(gr_x, gr_y, returnxy=True)
print '\nplotting...'
#mymap.contourf(grd[2], grd[3], DATA)
vmin = max(abs(np.min(DATA)), abs(np.max(DATA)))
if not ray_coverage:
    mymap.pcolormesh(grd[2], grd[3], DATA, cmap=tomo_colormap_2, vmin=-1*vmin/100., vmax=vmin/100.)
else:
    import matplotlib.cm as cm
    mymap.pcolormesh(grd[2], grd[3], DATA, cmap=cm.gray, vmax=10)
#plt.hexbin(grd[2], grd[3], DATA)
cbar = plt.colorbar(orientation='horizontal')
cbar.ax.tick_params(labelsize=12) 
plt.show()

# Clean 1.0% (dT)
import scipy.ndimage as ndimage
DATA_filt = ndimage.gaussian_filter(DATA, sigma=10.0, order=0)
#DATA_filt = DATA
mymap = Basemap(projection=projection, lon_0=long_0, lat_0=0)
mymap.drawcoastlines()
# this one is for CMB, based on your notes and corresponding to 1.%
mymap.pcolormesh(grd[2], grd[3], DATA_filt, cmap=tomo_colormap_2, vmin=-0.02537847, vmax=0.02537847)
cbar = plt.colorbar(orientation='horizontal')
cbar.ax.tick_params(labelsize=16)
cbar.ax.set_xticklabels(['-1.0%', ' ', ' ', ' ', '0%', ' ', ' ', ' ', '1.0%'])
plt.show()

# Clean 1.5%
import scipy.ndimage as ndimage
DATA_filt = ndimage.gaussian_filter(DATA, sigma=5.0, order=0)
DATA_filt = DATA
mymap = Basemap(projection=projection, lon_0=long_0, lat_0=0)
mymap.drawcoastlines()
mymap.pcolormesh(grd[2], grd[3], DATA_filt, cmap=tomo_colormap_2, vmin=-0.0380677, vmax=0.0380677)
#mymap.pcolormesh(grd[2], grd[3], DATA_filt, cmap=tomo_colormap_2, vmin=-0.02537848, vmax=0.02537848)
#plt.colorbar()
#plt.colorbar(orientation="horizontal")
cbar = plt.colorbar(orientation='horizontal')
cbar.ax.tick_params(labelsize=16)
cbar.ax.set_xticklabels(['-1.5%', ' ', ' ', ' ', '0%', ' ', ' ', ' ', '1.5%'])
plt.show()

# FOR RAY COVERAGE
import scipy.ndimage as ndimage
from matplotlib.colors import LogNorm
DATA_filt = ndimage.gaussian_filter(DATA, sigma=5.0, order=0)
DATA_filt = DATA
mymap = Basemap(projection=projection, lon_0=long_0, lat_0=0)
mymap.drawcoastlines(color='white')
vmin = max(abs(np.min(DATA)), abs(np.max(DATA)))
mymap.pcolormesh(grd[2], grd[3], DATA_filt, norm=LogNorm(vmin=0.1, vmax=330))
cbar = plt.colorbar(orientation='horizontal')
cbar.ax.tick_params(labelsize=16)
#cbar.ax.set_xticklabels(['0.1', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '285'])
plt.show()

# ==================== PICKLING THE REQUIRED INFORMATION ===============================
print "Pickling the final result..."
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
print "DONT FORGET TO MOVE THE %s TO ANOTHER DIRECTORY..." %(os.path.basename(processed_events_add) + '_MAP')

# For opening the files!
# 24 lines
#import glob
#from matplotlib.colors import LinearSegmentedColormap
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#from mpl_toolkits.basemap import Basemap
#import multiprocessing
#import numpy as np
#from obspy.core.util import locations2degrees
#import os
#import pickle
#import sys
#import time
#
#data_fi = open(os.path.join('.', 'DATA'), 'r')
#projec_fi = open(os.path.join('.', 'projection'), 'r')
#grd_fi = open(os.path.join('.', 'grd'), 'r')
#tomo_color_fi = open(os.path.join('.', 'tomo_colormap_2'), 'r')
#long_0_fi = open(os.path.join('.', 'long_0'), 'r')
#
#DATA = pickle.load(data_fi)
#projection = pickle.load(projec_fi)
#grd = pickle.load(grd_fi)
#tomo_colormap_2 = pickle.load(tomo_color_fi)
#long_0 = pickle.load(long_0_fi)
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
import scipy.ndimage as ndimage
from matplotlib.colors import LogNorm
DATA_filt = ndimage.gaussian_filter(DATA, sigma=5.0, order=0)
#DATA_filt = DATA
mymap = Basemap(projection=projection, lon_0=180, lat_0=0)
mymap.drawcoastlines(color='white')

#mymap.pcolormesh(grd[2], grd[3], DATA_filt, cmap=cm.gray, vmin=0, vmax=1.0)
mymap.pcolormesh(grd[2], grd[3], DATA_filt, norm=LogNorm(vmin=0.1, vmax=285))
#plt.colorbar()
plt.colorbar(orientation="horizontal")
cbar.ax.tick_params(labelsize=12) 
plt.show()

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
