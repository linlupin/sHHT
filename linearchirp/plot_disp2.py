####################################################################################
# To plot all the freq dispersion recorded in the hdf5 file (10000 trails recorded)
####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import ndimage

f = h5py.File("sHHTFD.hdf5", "r")
amp = f['/info/SNR'][()]                # SNR of the original data
input_noise_level = f['/info/INL'][()]  # Input noise level to perform sHHT
mdfreq_disp = f['/info/FD_sHHT'][()]    # Frequency dispersion between IFs obtained from sHHT (IF2) and conventional HHT (IF1)
mfreq_disp = f['/info/dFD_sHHT'][()]    # Frequency dispersion ratio between IF2 and IF1
nmdfreq_disp = f['/info/RFD_sHHT'][()]  # Frequency dispersion between IF obtained from sHHT (IF2) and real signal IF
nmfreq_disp = f['/info/RdFD_sHHT'][()]  # Frequency dispersion ratio between IF2 and real signal IF

label=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
fig=plt.figure(figsize=(20,10))
y, x = np.meshgrid(input_noise_level, amp)
Z = ndimage.gaussian_filter(mdfreq_disp, sigma=1.0, order=0)
CS = plt.contourf(x,y,Z)
plt.xscale('log')
plt.xlabel('SNR', fontsize=15)
plt.yticks(label)
plt.ylabel('Input Noise Level', fontsize=15)
plt.colorbar(CS)
plt.savefig('mddisp_contour.pdf',format="pdf", dpi=1200)
plt.close(fig)

fig=plt.figure(figsize=(20,10))
y, x = np.meshgrid(input_noise_level, amp)
CS = plt.contourf(x,y,mfreq_disp)
plt.xscale('log')
plt.xlabel('SNR', fontsize=15)
plt.yticks(label)
plt.ylabel('Input Noise Level', fontsize=15)
plt.colorbar(CS)
plt.savefig('mdisp_contour.pdf',format="pdf", dpi=1200)
plt.close(fig)

fig=plt.figure(figsize=(20,10))
y, x = np.meshgrid(input_noise_level, amp)
Z = ndimage.gaussian_filter(nmfreq_disp, sigma=1.0, order=0)
CS = plt.contourf(x,y,Z)
plt.xscale('log')
plt.xlabel('SNR', fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(label)
plt.yticks(fontsize=20)
plt.ylabel('Input Noise Level', fontsize=30)
plt.colorbar(CS)
plt.savefig('nmdisp_contour.eps',format="eps", dpi=1200)
plt.close(fig)


plt.figure(figsize=(20,10))
y, x = np.meshgrid(input_noise_level, amp)
Z = ndimage.gaussian_filter(nmdfreq_disp, sigma=0.8, order=0)
CS = plt.contourf(x,y,Z)
cb=plt.colorbar()
cb.ax.tick_params(labelsize=18)
plt.xscale('log')
xtick=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
xlabel=[0.1, 0.2, " ", " ", 0.5, " ", " ", " ", " ", 1.0, 2.0, " ", " ", 5.0, " ", " ", " ", " ", 10.0]
plt.xlabel('SNR', fontsize=30)
plt.xticks(xtick,labels=xlabel,fontsize=20)
plt.yticks(label)
plt.yticks(fontsize=20)
plt.ylabel('Input Noise Level', fontsize=30)
plt.savefig('nmddisp_contour.pdf',format="pdf", dpi=1200)
plt.close(fig)
#plt.show()



