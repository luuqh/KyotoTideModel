# KyoTide

A very high-resolution tide model developed at Kyoto University by Hung Q. Luu for the Tsugaru Strait. 

Its original source based on OGCM developed by Dr Yoichi Ishikawa. Barotropic tidal components were integrated into the model and used to simulate sea level change due to astronomical tide in the Tsugaru Strait which connects Sea of Japan with the Pacific Ocean. 

The model is written in FORTRAN and use Message Passing Interface (MPI) to run in parallel. Its default configuration on Kyoto University's Cray computers uses 4 nodes with 32 CPUs.
