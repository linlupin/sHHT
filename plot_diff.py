######################################################################################################
# To plot the frequency dispersion of a specific case recorded in the hdf5 file (1000 trials stacked)
#######################################################################################################

import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import ndimage

f = h5py.File("sin_disp.hdf5", "r")                 # To read the data with a stable signal
#f = h5py.File("linearchirp_disp.hdf5", "r")        # To read the data with a linear chirp signal
#f = h5py.File("expchirp_disp.hdf5", "r")           # To read the data with a exponential chirp signal
SNR = f['/info/SNR'][()]                            # Internal SNR of the original data
input_noise_level = f['/info/INL'][()]              # Input noise level to generate the sHHT
rfreq = f['/info/rfreq'][()]                        # IF of the original signal
OSF = f['/info/OSF'][()]                            # IF determined by the traditional HHT
SSF = f['/info/SSF'][()]                            # IF determined by sHHT
t = f['/info/time'][()]                             # Time
freq_center = f['/info/freq_center'][()]            # frequency
shht = f['/info/sHHT'][()]                          # methods to perform emd package

fig = plt.figure(figsize=(20,10))
plt.plot(t,rfreq,color='blue')
plt.plot(t,OSF,color='black')
plt.plot(t,SSF,color='red')
plt.xlabel('Time', fontsize=30)
plt.xticks(fontsize=20)
plt.ylabel('Frequency', fontsize=30)
plt.ylim(0,3)
plt.yticks(fontsize=20)
#plt.show()
plt.savefig('diff_HHT.pdf',format="pdf", dpi=1200)
#plt.savefig('diff_HHT_linearchirp.pdf',format="pdf", dpi=1200)
plt.close(fig)

levels = np.arange(0.0,0.04,0.008)
fig = plt.figure(figsize=(20,10))
Z = ndimage.gaussian_filter(shht, sigma=3.0, order=0)
plt.contour(t,freq_center,Z,levels,colors='g')
plt.plot(t,rfreq,color='blue')
plt.plot(t,SSF,color='red')
plt.imshow(shht,origin='lower',cmap='pink_r',extent=[np.min(t),np.max(t),np.min(freq_center),np.max(freq_center)], aspect='auto')
plt.clim(0,0.045)                               # boundary of the color bar needs to be revised with different data
plt.xlabel('Time', fontsize=30)
xtick=[0, 2, 4, 6, 8, 10]
xlabel=[0, 2, 4, 6, 8, 10]
plt.xticks(xtick,labels=xlabel,fontsize=20)
#ytick=[0, 1, 2, 3, 4, 5]
#ylabel=[0, 1, 2, 3, 4, 5]
ytick=[0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
ylabel=[0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
plt.ylabel('Frequency', fontsize=30)
plt.yticks(ytick,labels=ylabel,fontsize=20)
#plt.show()
plt.savefig('disp_sin.pdf',format="pdf", dpi=1200)
#plt.savefig('disp_linearchirp.pdf',format="pdf", dpi=1200)
plt.close(fig)

