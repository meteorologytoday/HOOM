#!/bin/bash


julia generate_SCRIP_format.jl \
    --input-file=/seley/tienyiah/CESM_domains/domain.ocn.gx3v7.120323.nc    \
    --output-file=domain.ocn.gx3v7.120323.SCRIP.nc    \
    --center-lon=xc     \
    --center-lat=yc     \
    --corner-lon=xv     \
    --corner-lat=yv

julia generate_SCRIP_format.jl \
    --input-file=/seley/tienyiah/CESM_domains/domain.lnd.fv4x5_gx3v7.091218.nc    \
    --output-file=domain.lnd.fv4x5_gx3v7.091218.SCRIP.nc    \
    --center-lon=xc     \
    --center-lat=yc     \
    --corner-lon=xv     \
    --corner-lat=yv     \
    --mask-flip

 mpirun -np 4 ESMF_RegridWeightGen -s domain.ocn.gx3v7.120323.SCRIP.nc -d domain.lnd.fv4x5_gx3v7.091218.SCRIP.nc -m conserve2nd -w w.nc --user_areas
