# Stacked Hilbert-Huang transform (sHHT)
This repository contains major python codes to perform the simulation in the following paper: "Effectiveness of Stacks in the Stacked Spectrum Hilbert-Huang Transform" - Yen al. (2024).

In the aforementioned paper, we descirbe the different frequency dispersions generated from conventional HHT and sHHT. We also demonstrate the frequency dispersions derived from 3 signal examples, including a stable sinudoidal wave, a linear chirp signal and a logarithmic chirp signal. At final, we also perform the sHHT to examine the known GW events, similar to the application of sHHT in [Hu et al. (2022)](https://arxiv.org/abs/2207.06714) to study the 11 GW (gravitational wave) events recorded before O3. In order to further confirm the practiability of the sHHT, we examined the sHHT on the GW events detected in O3 recorded in [GWTC-2](https://arxiv.org/abs/2010.14527) and [GWTC-3](https://arxiv.org/abs/2111.03606) (i.e., GW190814 and GW200129_065458). In comparison to both EMD (empirical mode decomposition) and ensemble empirical mode decomposition (EEMD) introduced on sHHT in [Hu et al. (2022)](https://arxiv.org/abs/2207.06714), here our codes are only based on EMD since [Hu et al. (2022)](https://arxiv.org/abs/2207.06714) has already demonstrated that a similar performance can be achieved by EMD or EEMD.

Our codes to perform the sHHT are based on the [EMD package 0.6.2](https://emd.readthedocs.io/en/stable/index.html) developed in python. EMD package is an open dource software basec on [Andrew J. Quinn (2021)](https://joss.theoj.org/papers/10.21105/joss.02977). To apply the EEMD, one can simply change the function from `emd.sift.sift` to `emd.sift.ensemble_sift`, and needs to inlcude a specific ensemble noise level (i.e., `ensemble_noise =` ).





