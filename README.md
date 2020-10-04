# KyoTide

A very high-resolution tide model developed at Kyoto University by Hung Q. Luu for the Tsugaru Strait. 

Its original source based on OGCM developed by Dr Yoichi Ishikawa. It was then incorporated barotropic tidal components from forcing NAO99Jb tidal output along open boundaries and Schwiderski (1980) astronomical parameters. It is used to simulate sea level change due to astronomical tide in the Tsugaru Strait which connects Sea of Japan with the Pacific Ocean. 

The model is written in FORTRAN and uses OpenMPI to run in parallel HPC. Its default configuration on Kyoto University's Cray computers uses 4 nodes with 32 CPUs.
