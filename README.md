# Stacked Hilbert-Huang transform (sHHT)
This repository contains major python codes to perform the simulation in the following paper: "Effectiveness of Stacks in the Stacked Spectrum Hilbert-Huang Transform" - Yen al. (2024).

In the aforementioned paper, we descirbe the different frequency dispersions generated from conventional HHT and sHHT. We also demonstrate the frequency dispersions derived from 3 signal examples, including a stable sinudoidal wave, a linear chirp signal and a logarithmic chirp signal. At final, we also perform the sHHT to examine the known GW events, similar to the application of sHHT in [Hu et al. (2022)](https://arxiv.org/abs/2207.06714) to study the 11 GW (gravitational wave) events recorded before O3. In order to further confirm the practiability of the sHHT, we examined the sHHT on the GW events detected in O3 recorded in [GWTC-2](https://arxiv.org/abs/2010.14527) and [GWTC-3](https://arxiv.org/abs/2111.03606) (i.e., GW190814 and GW200129_065458). In comparison to both EMD (empirical mode decomposition) and ensemble empirical mode decomposition (EEMD) introduced on sHHT in [Hu et al. (2022)](https://arxiv.org/abs/2207.06714), here our codes are only based on EMD since [Hu et al. (2022)](https://arxiv.org/abs/2207.06714) has already demonstrated the similar performance can be achieved by EMD and EEMD.



