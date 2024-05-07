#####################################################################################
#To show the stacked Hilbert spectrum of the GW data with different input noise level
#####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import h5py

matplotlib.use('TkAgg')

f = h5py.File("GW190814.hdf5", "r")     # We present GW190814 as one example
t = f['/info/TIME'][()]                 # Time series information  
gps = f['/info/gps'][()]                # mergering time of the event recorded on GWTC
sample_rate = f['/info/rate'][()]       # sampling rate of the time series
H1 = f['/info/H1_strain'][()]           # Hanford whitened strain data
L1 = f['/info/L1_strain'][()]           # Livingston whitened strain data

#To simply plot the data information
#plt.figure(figsize=(20,10))
#plt.plot(t,H1,color='red',linestyle='dashed')
#plt.plot(t,L1,color='blue')
#plt.xlabel('Time', fontsize=20)
#plt.ylabel('Whitened Strain', fontsize=20)
#plt.show()

print(gps)

import emd
import scipy
import math
from scipy import signal
from scipy import ndimage
method='hilbert'                # method to perform emd package
n=1000                          # 1000 trials to stack HHT
fmin=30.0                       # recommended frequency range to search for the signal
fmax=500
nfreq=1000                      # frequency grids
print(np.std(H1))

imf_opts = {'max_iters': 200, 'stop_method':'sd', 'sd_thresh': 0.005, 'energy_thresh': 500}
# We investigate different input noise level effected on the obtained Hilbert spectrum 
input_noise_level=[0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7,1,3,5,7]

for a in input_noise_level:
    input_noise=np.random.randn(len(t),n)*a*np.std(H1)
    shht=np.zeros((nfreq,len(t)))

    for i in np.arange(n):
        x=H1+input_noise[:,i]           # Hanford data
        imf = emd.sift.sift(x, imf_opts=imf_opts)
        IP, IF, IA = emd.spectra.frequency_transform(imf, sample_rate, method)
        freq_edges, freq_centres = emd.spectra.define_hist_bins(fmin, fmax, nfreq, 'linear')
        freq_center=(freq_edges[0:nfreq]+freq_edges[1:nfreq+1])*.5
        f, hht = emd.spectra.hilberthuang(IF[:, :], IA[:, :], freq_edges, mode='power', sum_time=False) # Use this line with all IMFs

        shht=shht+hht/n

    print('a = ',a,',median(shht) =',np.median(shht),',std(shht) = ',np.std(shht,ddof=1),',min(shht) =',np.min(shht))
    fig = plt.figure(figsize=(20,10))
    CS=plt.imshow(shht,origin='lower',cmap='pink_r',vmin=0.0,vmax=int(100.0*(np.median(shht)+2.5*np.std(shht,ddof=1)))/100.0,extent=[np.min(t),np.max(t),np.min(freq_center),np.max(freq_center)], aspect='auto')
    plt.colorbar(CS)
    plt.xlabel('Time', fontsize=30)
    plt.xlim(gps-0.35, gps+0.15)
    plt.xticks(fontsize=20)
    plt.ylabel('Frequency', fontsize=30)
    plt.yticks(fontsize=20)
    plt.savefig('GW190814_H1_'+str(a)+'.pdf',format="pdf", dpi=1200)
    plt.close(fig)


for a in input_noise_level:
    input_noise=np.random.randn(len(t),n)*a*np.std(L1)
    shht=np.zeros((nfreq,len(t)))

    for i in np.arange(n):
        x=L1+input_noise[:,i]            # Livingston data
        imf = emd.sift.sift(x, imf_opts=imf_opts)
        IP, IF, IA = emd.spectra.frequency_transform(imf, sample_rate, method)
        freq_edges, freq_centres = emd.spectra.define_hist_bins(fmin, fmax, nfreq, 'linear')
        freq_center=(freq_edges[0:nfreq]+freq_edges[1:nfreq+1])*.5
        f, hht = emd.spectra.hilberthuang(IF[:, :], IA[:, :], freq_edges, mode='power', sum_time=False) # Use this line with all IMFs

        shht=shht+hht/n

    print('a = ',a,',median(shht) =',np.median(shht),',std(shht) = ',np.std(shht,ddof=1))
    fig = plt.figure(figsize=(20,10))
    CS=plt.imshow(shht,origin='lower',cmap='pink_r',vmin=0.0,vmax=int(100.0*(np.median(shht)+2.5*np.std(shht,ddof=1)))/100.0,extent=[np.min(t),np.max(t),np.min(freq_center),np.max(freq_center)], aspect='auto')
    plt.colorbar(CS)
    plt.xlabel('Time', fontsize=30)
    plt.xlim(gps-0.35, gps+0.15)
    plt.xticks(fontsize=20)
    plt.ylabel('Frequency', fontsize=30)
    plt.yticks(fontsize=20)
    plt.savefig('GW190814_L1_'+str(a)+'.pdf',format="pdf", dpi=1200)
    plt.close(fig)
