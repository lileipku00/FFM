import commands
import numpy as np
import os
from obspy import read, Trace
from obspy.signal import cosTaper, util
from obspy.taup.taup import getTravelTimes
import scipy.io
import subprocess

try:
    from obspy.core.util import locations2degrees
except Exception, error:
    print '---------'
    print error
    print '---------'
    from obspy.taup.taup import locations2degrees


#--------------------------- convSTF_gauss --------------------------------
def convSTF_gauss(tr, sigma=30.):

    gauss = lambda (t, s): 1. / (2. * np.pi * s**2.)**.5 \
                           * np.exp(-1*(t**2)/(2*(s**2)))

    df = tr.stats.sampling_rate
    dt = 1./df

    t = np.linspace(0., sigma * 20., sigma * 20 * df + 1)
    stf = gauss((t-sigma*10, sigma))
    nstf = len(stf)

    nfft = util.nextpow2(max(nstf, tr.stats.npts)) * 2
    stff = np.fft.rfft(stf, n=nfft) * dt
    trf = np.fft.rfft(tr, n=nfft) * dt
    tr.data = np.fft.irfft(stff * trf)[sigma*10*df:sigma*10*df+len(tr.data)] * df

    return tr

#--------------------------- readSTF --------------------------------
def readSTF(add_stf):
    
    """
    Read the Source Time Function
    """
    
    #reading and creating the Source Time Function
    stf_matlab = scipy.io.loadmat(add_stf)
    stf = stf_matlab['STF']

    STF = np.array([stf[0][0]])
    for i in range(1, len(stf)):
        STF = np.append(STF, stf[i])
    
    stats = {'network': 'STF', 
             'station': 'STF', 
             'location': '',
             'channel': '01', 
             'npts': len(STF), 
             'sampling_rate': 10.0,
             'mseed' : {'dataquality': 'D'}}
    STF_tr = Trace(data=STF, header=stats)
    
    return STF_tr

#--------------------------- convSTF --------------------------------
def convSTF(tr, STF_tr):
    
    """
    STF*Trace
    """
    
    conv_tr = tr.copy()
    STF_tr.resample(tr.stats.sampling_rate)
    conv_tr.data = np.convolve(STF_tr.data, tr.data)
    
    return conv_tr
 
#--------------------------- preprocessing function -------------------------
def preprocessing(tr, usesac = False, tr_address=False, resample=False, 
                    rtr=False, rmean=False, taper=False, hfreq=False, lfreq=False):
    #XXX: remove the mean has not implemented
     
    if usesac:
        pwd = commands.getoutput('pwd')
        os.chdir(os.path.abspath(tr_address))
        tmp_name = 'TEMP-SAC-TRACE'
        tr.write(tmp_name, format = 'SAC')
        
        p = subprocess.Popen(['sac'],
                             stdout = subprocess.PIPE,
                             stdin  = subprocess.PIPE,
                             stderr = subprocess.STDOUT )
        s = 'read ' + tmp_name + '\n'
        if resample: s += 'interpolate delta ' + str(1./resample) + '\n'
        if rtr: s += 'rtrend\n'
        if rmean: s += 'rmean\n'
        if taper: s += 'taper' + '\n'
        if hfreq: s += 'lowpass BUTTER CORNER ' + str(hfreq) + '\n'
        if lfreq: s += 'highpass BUTTER CORNER ' + str(lfreq) + '\n'
        s += 'write ' + tmp_name + '\n' + 'quit\n'
        
        out = p.communicate(s)
        tr = read(tmp_name)[0]
        os.remove(tmp_name)
        os.chdir(pwd)

    else:
        if resample:
            tr.resample(resample)
        if rtr:
            tr_rtr = RTR(tr = tr, degree = 2)
            tr.data = tr_rtr
        if taper:
            tr.data *= cosTaper(len(tr.data), p=0.05)
        if hfreq:
            tr.filter('lowpass', freq=hfreq, corners=2)
        if lfreq:
            tr.filter('highpass', freq=lfreq, corners=2)
    return tr

#--------------------------- RTR --------------------------------
def RTR(tr, degree = 2):
    
    """
    Remove the trend by Fitting a linear function to the trace 
    with least squares and subtracting it
    """
    
    raw_f = tr.copy()

    t = []
    b0 = 0
    inc = []
    
    b = raw_f.stats['starttime']

    for i in range(0, raw_f.stats['npts']):
        inc.append(b0)
        b0 = b0+1.0/raw_f.stats['sampling_rate'] 
        b0 = round(b0, 4)
        
    A = np.vander(inc, degree)
    (coeffs, residuals, rank, sing_vals) = np.linalg.lstsq(A, raw_f.data)
    
    f = np.poly1d(coeffs)
    y_est = f(inc)
    rt_c = raw_f.data-y_est
    
    return rt_c

#----------------------------- time_window -----------------------------
class time_window:
    def __init__(self, tr):
        self.id='%s.%s.%s.%s' %(tr.stats.network, tr.stats.station,\
                                    tr.stats.location, tr.stats.channel)
        self.stats=tr.stats
    def epi_dist(self, req_phase='Pdiff', tb=10, ta=25, model='iasp91'):
        self.dist = locations2degrees(lat1 = self.stats.sac.evla, \
                        long1 = self.stats.sac.evlo, \
                        lat2 = self.stats.sac.stla, \
                        long2 = self.stats.sac.stlo)
        tt = getTravelTimes(self.dist, self.stats.sac.evdp, model=model)
        t_phase = -12345.0
        for tt_item in tt:
            if tt_item['phase_name'] == req_phase:
                t_phase = tt_item['time']
                phase_exist = 'Y'
                break
        if t_phase != -12345.0:
            t_before = t_phase - tb
            t_after = t_phase + ta
        else:
            phase_exist = 'N'
            t_before = t_phase
            t_after = t_phase
        return (phase_exist, t_phase, t_before, t_after, self.dist)


