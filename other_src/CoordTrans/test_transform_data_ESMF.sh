#!/bin/bash

julia transform_data_ESMF.jl  \
    --w-file=w.nc     \
    --s-file=LENS_piControl_oQ_3_f45_g37_NKOM_EKMAN_20.ocn.h.monthly.0021.nc  \
    --d-file=test.nc  \
    --vars=T,h_ML    \
    --x-dim=Nx        \
    --y-dim=Ny        \
    --z-dim=Nz_bone   \
    --t-dim=time  
