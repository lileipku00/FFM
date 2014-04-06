"""
This script is just a 'TRASH' type script!
i.e.: manually change the script and run it for different cases
Goal: plot different meanall_dts that have been obtained from dispersion_dT_dA.py
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import pickle

parent_add = '/import/neptun-radler/hosseini-downloads/KASRA/SCRIPTS/gitHUB/myrepo_gitHUB/FFM/STATISTICS' + \
             '/statistics/FFM_PAPER_April_2014'

case_1_fio = open(os.path.join(os.path.join(parent_add, 'meanall_dt_GSN_1-5')))
case_1 = pickle.load(case_1_fio)


plt.ylabel('Time difference (dT)', fontsize='x-large', weight='bold')
plt.xlabel('Dominant Period', fontsize='x-large', weight='bold')
x = [2.7, 3.7, 5.3, 7.5, 10.6, 15.0, 21.2, 30.0]
plt.vlines(x, -10.0, 10.0)
plt.xlim(xmin=0.0)
plt.ylim(-1.0, 1.0)
plt.xticks(x, fontsize='x-large', weight='bold')
plt.yticks(np.arange(-1.0, 1.0, 0.05), fontsize='x-large', weight='bold')

plt.plot(x[::-1], case_1, 'b', linewidth=3, label='Pdiff, NLAMBDA=1.5 (GSN)')

plt.legend(prop={'size': 22})
plt.show()

