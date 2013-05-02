#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  all_plot_dT_dA.py
#   Purpose:   run plot_dT_dA.py for all available events
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import sys

# ------------------- INPUT -----------------------------
plt_all = False
# -------------------------------------------------------

bands = sys.argv[1]
remote_dir = '/import/neptun-radler/hosseini-downloads/KASRA/FFM' 

proc_ev_ls = glob.glob(os.path.join(remote_dir, '*.*.*.*'))
print '%s processed events found!' %(len(proc_ev_ls))
failed = 0
if os.path.isfile(os.path.join('.', 'mean_values.txt')):
    print 'mean_values.txt found in the directory...'
else:
    print '\nERRORS:'
    for i in range(len(proc_ev_ls)):
        print '.',
        sys.stdout.flush()
        try:
            evadd = proc_ev_ls[i]
            subprocess.check_call(['python', 'plot_dT_dA.py', bands, evadd])
        except Exception, e:
            print e
            failed += 1

print '\nplotting!'
fio = open(os.path.join('.', 'mean_values.txt'))
fi = fio.readlines()
x_all = []
y_all = []
weight_all = []
for i in range(len(fi)):
    fi[i] = fi[i].split(',')
    if np.isnan(float(fi[i][1])) or np.isnan(float(fi[i][2])):
        continue
    x = (float(fi[i][5]), float(fi[i][6]))
    y = (float(fi[i][1])*float(fi[i][5])+float(fi[i][2]),
          float(fi[i][1])*float(fi[i][6])+float(fi[i][2]))
    x_all.append(x)
    y_all.append(y)
    weight_all.append(float(fi[i][7]))
    if plt_all:
        plt.plot(x, y)
y1_all = np.array([])
y2_all = np.array([])
for i in range(len(y_all)):
    y1_all = np.append(y1_all, y_all[i][0]*weight_all[i]) 
    y2_all = np.append(y2_all, y_all[i][1]*weight_all[i])

plt.plot(x, (np.sum(y1_all)/sum(weight_all), 
            np.sum(y2_all)/sum(weight_all)), '--',
            linewidth = 3)

plt.ylabel('Time difference (dT)', fontsize = 'x-large', weight = 'bold')
plt.xlabel('Dominant Period', fontsize = 'x-large', weight = 'bold')
plt.xticks(fontsize = 'x-large', weight = 'bold')
plt.yticks(fontsize = 'x-large', weight = 'bold')
pltitle = 'All the processed events'
pltitle += '\nbands: %s' %(bands)
plt.title(pltitle, fontsize = 'x-large', weight = 'bold')

plt.show()
