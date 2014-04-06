'''
This script is just a 'trash' type script!
i.e.: manually change the script and run it for different cases
Goal: plot different meanall_dt that obtained from plot_dT_dA.py
'''
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle

case_1_fio = open(os.path.join('PRE_INVERSION_PAPER', 'CHANGE_BM', 'meanall_dt_selected_events'))
case_1 = pickle.load(case_1_fio)
case_2_fio = open(os.path.join('PRE_INVERSION_PAPER', 'meanall', 'meanall_dt_all'))
case_2 = pickle.load(case_2_fio)
case_3_fio = open(os.path.join('PRE_INVERSION_PAPER', 'CHANGE_BM', 'meanall_dt_selected_new_BM'))
case_3 = pickle.load(case_3_fio)
case_4_fio = open(os.path.join('PRE_INVERSION_PAPER', 'meanall', 'meanall_dt_selected_stations'))
case_4 = pickle.load(case_4_fio)
case_5_fio = open(os.path.join('PRE_INVERSION_PAPER', 'CHANGE_BM', 'meanall_dt_selected_30_BM'))
case_5 = pickle.load(case_5_fio)
case_6_fio = open(os.path.join('PRE_INVERSION_PAPER', 'NLAMBDA_1-5', 'meanall_dt_normal'))
case_6 = pickle.load(case_6_fio)
case_7_fio = open(os.path.join('PRE_INVERSION_PAPER', 'NLAMBDA_1-5', 'meanall_dt_10'))
case_7 = pickle.load(case_7_fio)
case_8_fio = open(os.path.join('PRE_INVERSION_PAPER', 'NLAMBDA_1-5', 'meanall_dt_30'))
case_8 = pickle.load(case_8_fio)

plt.ylabel('Time difference (dT)', fontsize = 'x-large', weight = 'bold')
plt.xlabel('Dominant Period', fontsize = 'x-large', weight = 'bold')
x = [2.7, 3.7, 5.3, 7.5, 10.6, 15.0, 21.2, 30.0]
plt.vlines(x, 0.0, 1.0)
plt.xlim(xmin=0.0)
plt.ylim(0.25, 0.8)
plt.ylim(0., 0.8)
plt.xticks(x, fontsize = 'x-large', weight = 'bold')
plt.yticks(np.arange(0.0, 0.9, 0.05), fontsize = 'x-large', weight = 'bold')

plt.plot(x[::-1], case_2, 'k', linewidth=3, label='All events')
plt.plot(x[::-1], case_4, 'r', linewidth=3, label='Selected stations')
plt.plot(x[::-1], case_1, 'b', linewidth=3, label='Selected events')
plt.plot(x[::-1], case_3, 'b--', linewidth=3, label='Selected events (modified 10%)')
plt.plot(x[::-1], case_5, 'b', linestyle='dashdot', linewidth=3, label='Selected events (modified 30%)')
plt.plot(x[::-1], case_6, 'g', linewidth=3, label='Selected events (NLAMBDA=1.5)')
plt.plot(x[::-1], case_7, 'g--', linewidth=3, label='Selected events (modified 10%, NLAMBDA=1.5)')
plt.plot(x[::-1], case_8, 'g', linestyle='dashdot', linewidth=3, label='Selected events (modified 30%, NLAMBDA=1.5)')
plt.legend(prop={'size':22})
plt.show()

