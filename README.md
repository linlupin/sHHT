# Stacked Hilbert-Huang transform (sHHT)
This repository contains major python codes to perform the simulation in the following paper: "Effectiveness of Stacks in the Stacked Spectrum Hilbert-Huang Transform" - Yen al. (2024).

In the aforementioned paper, we descirbe the different frequency dispersions generated from conventional HHT and sHHT. We also demonstrate the frequency dispersions derived from 3 signal examples, including a stable sinudoidal wave, a linear chirp signal and an exponential chirp signal. At final, we also perform the sHHT to examine the known GW events, similar to the application of sHHT in [Hu et al. (2022)](https://arxiv.org/abs/2207.06714) to study the 11 GW (gravitational wave) events recorded before O3. In order to further confirm the practiability of the sHHT, we examined the sHHT on the GW events detected in O3 recorded in [GWTC-2](https://arxiv.org/abs/2010.14527) and [GWTC-3](https://arxiv.org/abs/2111.03606) (i.e., GW190814 and GW200129_065458). In comparison to both EMD (empirical mode decomposition) and ensemble empirical mode decomposition (EEMD) introduced on sHHT in [Hu et al. (2022)](https://arxiv.org/abs/2207.06714), here our codes are only based on EMD since [Hu et al. (2022)](https://arxiv.org/abs/2207.06714) has already demonstrated that a similar performance can be achieved by EMD or EEMD.

Our codes to perform the sHHT are based on the [EMD package 0.6.2](https://emd.readthedocs.io/en/stable/index.html) developed in python. EMD package is an open dource software basec on [Andrew J. Quinn (2021)](https://joss.theoj.org/papers/10.21105/joss.02977). To apply the EEMD, one can simply change the function from `emd.sift.sift` to `emd.sift.ensemble_sift`, and needs to inlcude a specific ensemble noise level (i.e., `ensemble_noise =` ).
Here we can compare the conventional HHT with EEMD and sHHT using the following two flow charts.
![shht](https://github.com/linlupin/sHHT/blob/main/sHHT.png) 

In our investigations, we considered a stationary signal @ 1 Hz and chirp signals in linear and exponential behavior. We assume a sinusoidal function of sin(2πft) while f=0.1*(t+0.01) to simulate a linear chirp signal and f=2^{0.1*t} to simulate an exponential chirp signal. The `stationary`, `linearchirp` and `expchirp` folders are used to denote codes to inverstigate the signal of corresponding behaviors, respectively. Each folder has following 3 codes:  
-  ***Freq_disp.py***: To investigate the frequency disperison variation between the frequency of the real signal and the instantaneous frequency (IF) determined by sHHT following the change of different input noise levels and internal data SNRs.  
-  ***plot_disp2.py***: To plot the data stored by Freq_disp.py.  
-  ***diff_sHHT_xxx.py***: To examine the relation of the real signal frequency, IF determined by conventional HHT (IF1) and sHHT (IF2). `xxx` denotes the behavior of the signal.

***plot_diff.py*** in the root folder is used to plot the result stored by diff_sHHT_xxx.py.  
The GW folder includes the codes to deal with the data reduction of the gravitational wave (GW) data and to present the GW signal on the stacked Hilbert spectra. There are 4 files including in this folder:  
+  ***GW190814.ipynb***: Google Jupyter notebook to deal with the reduction of GW190814 strain data. The Q-transform plot is also included.  
+  ***GW200129.ipynb***: Google Jupyter notebook to deal with the reduction of GW200129_065458 strain data. The Q-transform plot is also included.  
More details about the data reduction of GW events can be referred to the tutorial resource provided by Gravitational Wave Open Science Center ([GWOSC](https://colab.research.google.com/github/gw-odw/)). (One can choose the tutorial resources from the latest ODW to practice)
+ ***GW_sHHT_lv.py***: To present the stacked Hilbert spectra of a specific GW event generated with different input noiselevels. One can fix a specific level to gain the best visual plot to generate the stacked Hilbert spectrum using the following code.     
+ ***GW_sHHT.py***: To present the stacked Hilbert spectra of a specific GW event generated with a specific input noise level.  
