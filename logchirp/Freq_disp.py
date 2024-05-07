# -*- coding: utf-8 -*-
##################################################################
# To count for the frequency dispersion of the linear chirp signal
##################################################################

import emd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
import h5py
import scipy
import math
from scipy import signal

matplotlib.use('Agg') 

ff=h5py.File("sHHTFD.hdf5","w")

data_noise_level=10.0    # internal noise level
#method='nht'            
method='hilbert'         # method to perform emd
n=500                    # number of HHT to stack (sHHT)
nfreq=1000               # number of frequency grids
# Frequency range
fmin=0.01
fmax=5.0
ind = 0                  # index to count for the simulation loop

t=np.arange(0,10,0.01)       # Time from 0 - 10 sec
freq=pow(2,0.1*t)            # 2^{0.1*t}; Real frequency of the signal is d(freq*t)/dt   
rfreq=pow(2,0.1*t)*(1+0.1*t*math.log(2))   # d(freq*t)/dt ; IF of the original signal
omega=2.0*np.pi*freq         # Angular Frequency
# SNR of one data can be treated as "amp/data_noise_level"
amp=np.arange(1.0,101.0,1.0) # Amplitude level(range) of the signal
# Additional input noise level injected on the data to perform the sHHT
input_noise_level=np.arange(0.01,1.01,0.01)  

noise=np.random.randn(len(t))
sample_rate=1/(t[1]-t[0])     # Sampling rate

OSF = np.zeros( len(t) )      # instantaneous frequency obtained from conventional HHT (IF1)
SSF = np.zeros( len(t) )      # instantaneous frequency obtained from sHHT (IF2)
smooth = np.zeros( len(t) )   # To smooth the obtained IF
Ofreq_disp =  np.zeros( len(t) )   # Frequency dispersion between the IF1 and real signal IF
freq_disp =  np.zeros( len(t) )    # Frequency dispersion ratio between IF2 and IF1
dfreq_disp =  np.zeros( len(t) )   # Frequency dispersion between the IF2 and IF1
nfreq_disp =  np.zeros( len(t) )   # Frequency dispersion ratio between IF2 and real signal IF
ndfreq_disp =  np.zeros( len(t) )  # Frequency dispersion between the IF2 and real signal IF
# standard deviation of the freq_disp
sfreq_disp =  np.zeros( (len(amp), len(input_noise_level)) )
# mean value of the freq_disp
mfreq_disp =  np.zeros( (len(amp), len(input_noise_level)) )
# mean value of the dfreq_disp
mdfreq_disp =  np.zeros( (len(amp), len(input_noise_level)) )
# mean value of the nfreq_disp
nmfreq_disp =  np.zeros( (len(amp), len(input_noise_level)) )
# mean value of the ndfreq_disp
nmdfreq_disp =  np.zeros( (len(amp), len(input_noise_level)) )



max_freqpow = None
freq_idx = None

start_time = time.time()
print('Start!')


imf_opts = {'max_iters': 200, 'stop_method':'sd', 'sd_thresh': 0.005, 'energy_thresh': 50}
max = 0

for i in np.arange(len(amp)):
    x=np.sin(omega*t)+noise*data_noise_level/amp[i]
    imf = emd.sift.sift(x, imf_opts=imf_opts)
    OIP, OIF, OIA = emd.spectra.frequency_transform(imf, sample_rate, method)
    freq_edges, freq_centres = emd.spectra.define_hist_bins(fmin, fmax, nfreq, 'linear')
    freq_center=(freq_edges[0:nfreq]+freq_edges[1:nfreq+1])*.5
    f, hht = emd.spectra.hilberthuang(OIF[:, :], OIA[:, :], freq_edges, mode='power', sum_time=False) # Use this line with all IMFs

    for ii in np.arange(len(t)):
         a = hht[:,ii]
         for idx, num in enumerate(a):
             if max_freqpow is None or num > max_freqpow:
                   max_freqpow = num
                   freq_idx = idx
         OSF[ii]=freq_center[freq_idx]
#         print(freq_idx)
         max_freqpow=None
    smooth=signal.savgol_filter(OSF,100,2)
    for ii in np.arange(len(t)-1):
         # To smooth the singularity when a weird peak occurs
         if (OSF[ii+1]-OSF[ii]) > 0.5 or (OSF[ii+1]-OSF[ii]) < -0.5:
             OSF[ii]=smooth[ii]
             OSF[ii+1]=smooth[ii+1]
         Ofreq_disp[ii]=OSF[ii]-rfreq[ii]


    for w in np.arange(len(input_noise_level)):
        shht=np.zeros((nfreq,len(t)))
        for times in np.arange(n):
            input_noise=np.random.randn(len(t),n)*input_noise_level[w]
            x1=x+input_noise[:,times]
            imf = emd.sift.sift(x1, imf_opts=imf_opts)
            IP, IF, IA = emd.spectra.frequency_transform(imf, sample_rate, method)
            freq_edges, freq_centres = emd.spectra.define_hist_bins(fmin, fmax, nfreq, 'linear')
            freq_center=(freq_edges[0:nfreq]+freq_edges[1:nfreq+1])*.5
            f, hht = emd.spectra.hilberthuang(IF[:, :], IA[:, :], freq_edges, mode='power', sum_time=False) # Use this line with all IMFs
#####
            shht=shht+hht/n

        for ii in np.arange(len(t)):
            a = shht[:,ii]
            for idx, num in enumerate(a):
                if max_freqpow is None or num > max_freqpow:
                     max_freqpow = num
                     freq_idx = idx
            SSF[ii]=freq_center[freq_idx]
            max_freqpow=None
        smooth=signal.savgol_filter(SSF,100,2)

        for ii in np.arange(len(t)-1):
              # To smooth the singularity when a weird peak occurs
             if (SSF[ii+1]-SSF[ii]) > 0.5 or (SSF[ii+1]-SSF[ii]) < -0.5:
                 SSF[ii]=smooth[ii]
                 SSF[ii+1]=smooth[ii+1]
             dfreq_disp[ii]=SSF[ii]-OSF[ii]
             ndfreq_disp[ii]=SSF[ii]- rfreq[ii]
             freq_disp[ii]=abs(SSF[ii]-OSF[ii])/OSF[ii]
             nfreq_disp[ii]=abs(SSF[ii]-rfreq[ii])/rfreq[ii]

        # To count for the average and standard deviation of the freq_disp at a specific time
        # We ignore the results obtained in the boundaries (< 1 s and > 9 s)
        mfreq_disp[i,w]= np.mean(freq_disp[100:900])   
        mdfreq_disp[i,w]= np.mean(abs(dfreq_disp[100:900]))   
        nmfreq_disp[i,w]= np.mean(nfreq_disp[100:900])   
        nmdfreq_disp[i,w]= np.mean(abs(ndfreq_disp[100:900]))  
        sfreq_disp[i,w]= np.std(freq_disp[100:900])
       
         # To check the first result to yield large frequency dispersion 
        if ind==0 and mfreq_disp[i,w] > 1.0:
              fig=plt.figure(figsize=(20,10))
              plt.plot(t[10:990],rfreq[10:990],color='black')
              plt.plot(t[10:990],OSF[10:990],color='red')
              plt.plot(t[10:990],SSF[10:990],color='blue')
              plt.xlabel('Time', fontsize=15)
              plt.ylabel('Frequency', fontsize=15)
              plt.savefig('f_t.pdf',format="pdf", dpi=1200)
              plt.close(fig)
    
        if ind==0 and mfreq_disp[i,w] > 1.0:
             fig=plt.figure(figsize=(20,10))
             plt.plot(t[10:990],Ofreq_disp[10:990],color='red')
             plt.plot(t[10:990],dfreq_disp[10:990],color='blue')
             plt.xlabel('Time', fontsize=15)
             plt.ylabel('Frequency dispersion', fontsize=15)
             plt.savefig('disp.pdf',format="pdf", dpi=1200)
             plt.close(fig)
             print(amp[i],input_noise_level[w],mfreq_disp[i,w])
             ind = ind + 1


label=[1,2,3,4,5,6,7,8,9,10]
print(amp.shape)

fig=plt.figure(figsize=(20,10))
y, x = np.meshgrid(input_noise_level, amp/data_noise_level)
CS = plt.contourf(x,y,mdfreq_disp)
plt.xscale('log')
plt.xlabel('SNR', fontsize=15)
plt.yticks(label)
plt.ylabel('Input Noise Level', fontsize=15)
plt.colorbar(CS)
plt.savefig('mddisp_contour.pdf',format="pdf", dpi=1200)
plt.close(fig)
fig=plt.figure(figsize=(20,10))

y, x = np.meshgrid(input_noise_level, amp/data_noise_level)
CS = plt.contourf(x,y,mfreq_disp)
plt.xscale('log')
plt.xlabel('SNR', fontsize=15)
plt.yticks(label)
plt.ylabel('Input Noise Level', fontsize=15)
plt.colorbar(CS)
plt.savefig('mdisp_contour.pdf',format="pdf", dpi=1200)
plt.close(fig)

print(amp.shape)
fig=plt.figure(figsize=(20,10))
y, x = np.meshgrid(input_noise_level, amp/data_noise_level)
CS = plt.contourf(x,y,nmdfreq_disp)
plt.xscale('log')
plt.xlabel('SNR', fontsize=15)
plt.yticks(label)
plt.ylabel('Input Noise Level', fontsize=15)
plt.colorbar(CS)
plt.savefig('nmddisp_contour.pdf',format="pdf", dpi=1200)
plt.close(fig)

fig=plt.figure(figsize=(20,10))
y, x = np.meshgrid(input_noise_level, amp/data_noise_level)
CS = plt.contourf(x,y,nmfreq_disp)
plt.xscale('log')
plt.xlabel('SNR', fontsize=15)
plt.yticks(label)
plt.ylabel('Input Noise Level', fontsize=15)
plt.colorbar(CS)
plt.savefig('nmdisp_contour.pdf',format="pdf", dpi=1200)
plt.close(fig)

print("information!")
print(f.ndim)
print(f.shape)

# Output the data of the results
f1=ff.create_group("info")
f1["SNR"] = amp/data_noise_level   # original SNR of the data
f1["INL"] = input_noise_level      # input noise level to perform sHHT
f1["FD_sHHT"] = mdfreq_disp        # Frequency dispersion between IF2 and IF1
f1["dFD_sHHT"] = mfreq_disp        # Frequency dispersion ratio between IF2 and IF1
f1["RFD_sHHT"] = nmdfreq_disp      # Frequency dispersion between the IF2 and real IF of the signal
f1["RdFD_sHHT"] = nmfreq_disp      # Frequency dispersion ratio between IF2 and real IF of the signal

ff.close()



end_time = time.time()
print("Elapsed Time: %f s" % (end_time - start_time))

