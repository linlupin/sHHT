#####################################################
#To show the stacked Hilbert spectrum of the GW data 
#####################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import h5py

matplotlib.use('TkAgg')

event = 'GW190814'                    # event name to label the file name
extname = '.hdf5'                     # extensive file name
event_file = event+extname            # file name

f = h5py.File(event_file, "r")
t = f['/info/TIME'][()]               # Time series information
gps = f['/info/gps'][()]              # mergering time of the event recorded on GWTC
sample_rate = f['/info/rate'][()]     # sampling rate of the time series
H1 = f['/info/H1_strain'][()]         # Hanford whitened strain data 
L1 = f['/info/L1_strain'][()]         # Livingston whitened strain data

#To simply plot the data information
plt.figure(figsize=(20,10))
plt.plot(t,H1,color='red',linestyle='dashed')
plt.plot(t,L1,color='blue')
plt.xlabel('Time', fontsize=20)
plt.ylabel('Whitened Strain', fontsize=20)
plt.show()

print(gps)

import emd
import scipy
import math
from scipy import signal
from scipy import ndimage
method='hilbert'                      # method to perform emd package
n=10000                               # 10000 trials to stack HHT
input_noise_level=0.9                 # Best/Specific (?) input noise level to geenrate the Hilbert spectrum
fmin=30.0                             # recommended frequency range to search for the signal
fmax=1000.0
nfreq=1000                            # frequency grids
print(np.std(L1))

imf_opts = {'max_iters': 200, 'stop_method':'sd', 'sd_thresh': 0.005, 'energy_thresh': 500}
shht=np.zeros((nfreq,len(t)))
input_noise=np.random.randn(len(t),n)*input_noise_level*np.std(L1)
freq_center = np.zeros(nfreq)      
freq_edges = np.zeros(nfreq+1)


for i in np.arange(n):
    x=L1+input_noise[:,i]              # Livingston data; for Hanford data, please change to H1
    imf = emd.sift.sift(x, imf_opts=imf_opts)
    IP, IF, IA = emd.spectra.frequency_transform(imf, sample_rate, method)
    freq_edges, freq_centres = emd.spectra.define_hist_bins(fmin, fmax, nfreq, 'linear')
    freq_center=(freq_edges[0:nfreq]+freq_edges[1:nfreq+1])*.5
    f, hht = emd.spectra.hilberthuang(IF[:, :], IA[:, :], freq_edges, mode='power', sum_time=False) # Use this line with all IMFs
    shht=shht+hht/n


print(np.median(shht))                 # To check the possible distribution of the Hilbert power/energy 
print(np.std(shht))
#levels = np.arange(0.0,0.1,0.025)
fig = plt.figure(figsize=(20,10))
#Z = ndimage.gaussian_filter(shht, sigma=2.0, order=0)
#plt.contour(t,freq_center,shht,levels,colors='g')
#CS=plt.contour(t[5000:np.size(H1)-5000],freq_center,shht,levels,cmap='pink_r')
#plt.colorbar(CS)
gps_i = int(gps/10.0)                  # Integer part of the merging time
gps_10i = gps_i * 10
plt.imshow(shht,origin='lower',cmap='pink_r',extent=[np.min(t),np.max(t),np.min(freq_center),np.max(freq_center)], aspect='auto')
#plt.clim([0.0,int(100.0*(np.median(shht)+10.0*np.std(shht,ddof=1)))/100.0])
plt.clim([0.0,0.020])  #GW190814
cb=plt.colorbar()
cb.ax.tick_params(labelsize=18) 
#plt.text(gps-0.095,430,'(a)',{'fontsize':40})
plt.text(gps-0.325,430,'(b)',{'fontsize':40})
#plt.contourf(t,freq_center,Z,cmap='pink_r')
gps_i = int(gps/10.0)
plt.xlabel('Time (+'+ str(gps_10i) + ' sec)' , fontsize=30)
plt.xlim(gps-0.35, gps+0.15)
#plt.xlim(gps-0.1, gps+0.1)
plt.xticks(fontsize=20)
plt.ylabel('Frequency', fontsize=30)
#plt.yscale("log")
xtick=[gps-0.3, gps-0.2, gps-0.1, gps, gps+0.1]
xlabel=[gps-gps_10i-0.3, gps-gps_10i-0.2,gps-gps_10i-0.1, gps-gps_10i, gps-gps_10i+0.1]
#xtick=[gps-0.1, gps-0.075, gps-0.05, gps-0.025, gps, gps+0.025, gps+0.05, gps+0.075, gps+0.1]
#xlabel=[round(gps-gps_10i-0.1,3), round(gps-gps_10i-0.075,3), round(gps-gps_10i-0.05,3), round(gps-gps_10i-0.025,3), round(gps-gps_10i,3), round(gps-gps_10i+0.025,3), round(gps-gps_10i+0.05,3), round(gps-gps_10i+0.075,3), round(gps-gps_10i+0.1,3)]
plt.xticks(xtick,labels=xlabel,fontsize=20)
plt.xticks(fontsize=20)
ytick=[100, 200, 300, 400, 500]
ylabel=[100, 200, 300, 400, 500]
plt.yticks(ytick,labels=ylabel,fontsize=20)
plt.ylim(30, 500)
plt.yticks(fontsize=20)
plt.savefig(event+'_L1_'+str(input_noise_level)+'.pdf',format="pdf", dpi=1200)
plt.close(fig)
#plt.show()
