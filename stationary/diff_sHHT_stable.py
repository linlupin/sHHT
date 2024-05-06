###################################################################
# To plot for the IF of traditional HHT and sHHT for a stable signal
###################################################################

import emd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
import h5py
import scipy
import math
from scipy import signal
from scipy import ndimage

matplotlib.use('TkAgg')
ff=h5py.File("sin_disp2.hdf5","w")    # output the result for a plot

data_noise_level=10.0                 # internal noise level of the data
amp = 20.0                            # amplitude of the signal
                                      # The orignal data SNR can be treated as amp/data_noise_level
#method='nht'
method='hilbert'                      # methods to perform emd package
n=1000                                # stack HHT for 1000 times
#input_noise_level=0.3
input_noise_level=0.01                # input noise level to perform sHHT 


nfreq=600                             # frequency grids
fmin=0.01                             # frequency range: 0.01 - 3 Hz
fmax=3.0

t=np.arange(0,10,0.01)
freq=1.0                              # stable srequency @ 1Hz
rfreq=np.ones(len(t))                 # d(freq*t)/dt (IF of the original signal)
omega=2.0*np.pi*freq                  # Angluar frequency

noise=np.random.randn(len(t))*data_noise_level/amp    # random noise         
sample_rate=1/(t[1]-t[0])             # sampling rate  

start_time = time.time()
print('Start!')

imf_opts = {'max_iters': 200, 'stop_method':'sd', 'sd_thresh': 0.005, 'energy_thresh': 50}
max = 0
shht=np.zeros((nfreq,len(t)))
print(input_noise_level)
input_noise=np.random.randn(len(t),n)*input_noise_level

max_freqpow = None
freq_idx = None
OSF = np.zeros( len(t) )              # instantaneous frquency (IF) determined by HHT (before stacked)

for i in np.arange(n):
    x=np.sin(omega*t)+noise+input_noise[:,i]
    imf = emd.sift.sift(x, imf_opts=imf_opts)
    OIP, OIF, OIA = emd.spectra.frequency_transform(imf, sample_rate, method)
    freq_edges, freq_centres = emd.spectra.define_hist_bins(fmin, fmax, nfreq, 'linear')
    freq_center=(freq_edges[0:nfreq]+freq_edges[1:nfreq+1])*.5
    f, hht = emd.spectra.hilberthuang(OIF[:, :], OIA[:, :], freq_edges, mode='power', sum_time=False) # Use this line with all IMFs

    shht=shht+hht/n
   
    if i==0:         # The first trial is used to determined the result obtained from the traditional HHT   
        for ii in np.arange(len(t)):
             a = hht[:,ii]
             for idx, num in enumerate(a):
                  if max_freqpow is None or num > max_freqpow:
                     max_freqpow = num
                     freq_idx = idx
             OSF[ii]=freq_center[freq_idx]
             max_freqpow=None
        smooth=signal.savgol_filter(OSF,100,2)

        # To smooth the singularity when a weird peak occurs
        for ii in np.arange(len(t)-1):
             if (OSF[ii+1]-OSF[ii]) > 0.5 or (OSF[ii+1]-OSF[ii]) < -0.5:
                 OSF[ii]=smooth[ii]
                 OSF[ii+1]=smooth[ii+1]


max_freqpow = None
freq_idx = None
SSF = np.zeros( len(t) )              # IF determined by sHHT

for ii in np.arange(len(t)):
    a = shht[:,ii]
    for idx, num in enumerate(a):
        if max_freqpow is None or num > max_freqpow:
             max_freqpow = num
             freq_idx = idx
    SSF[ii]=freq_center[freq_idx]
#   print(freq_idx)
    max_freqpow=None
smooth=signal.savgol_filter(SSF,100,2)

# To smooth the singularity when a weird peak occurs
for ii in np.arange(len(t)-1):
    if (SSF[ii+1]-SSF[ii]) > 0.5 or (SSF[ii+1]-SSF[ii]) < -0.5:
        SSF[ii]=smooth[ii]
        SSF[ii+1]=smooth[ii+1]

# To save the data for further plot and confirmation
f1=ff.create_group("info")
f1["SNR"] = data_noise_level/amp      # Internal SNR of the original data
f1["INL"] = input_noise_level         # Input noise level to generate the sHHT
f1["rfreq"] = rfreq                   # IF of the original signal
f1["OSF"] = OSF                       # IF determined by the traditional HHT
f1["SSF"] = SSF                       # IF determined by sHHT
f1["time"] = t                        # Time
f1["freq_center"] = freq_center       # frequency
f1["sHHT"] = shht                     # methods to perform emd package
ff.close()

# To generate the plot of the IF of the signal, IF determined by conventional HHT and sHHT
fig = plt.figure(figsize=(20,10))
plt.plot(t,rfreq,color='blue')
plt.plot(t,OSF,color='black')
plt.plot(t,SSF,color='red')
plt.xlabel('Time', fontsize=30)
plt.xticks(fontsize=20)
plt.ylabel('Frequency', fontsize=30)
plt.ylim(0,3)
plt.yticks(fontsize=20)
#plt.savefig('diff_HHT2.pdf',format="pdf", dpi=1200)
#plt.close(fig)
plt.show()

# To generate a plot of IF determined by conventional HHT and sHHT imposed on the Hilbert spectrum
levels = np.arange(0.0,0.1,0.025)
fig = plt.figure(figsize=(20,10))
Z = ndimage.gaussian_filter(shht, sigma=3.0, order=0)
plt.contour(t,freq_center,Z,levels,colors='g')
plt.plot(t,rfreq,color='blue')
plt.plot(t,SSF,color='red')
CS=plt.contourf(t,freq_center,shht,cmap='pink_r')
plt.colorbar(CS)
plt.xlabel('Time', fontsize=30)
plt.xticks(fontsize=20)
plt.ylabel('Frequency', fontsize=30)
plt.yticks(fontsize=20)
#plt.savefig('disp_sin2.pdf',format="pdf", dpi=1200)
#plt.close(fig)
plt.show()
