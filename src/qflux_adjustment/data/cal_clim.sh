#!/bin/bash

d=/seley/tienyiah/cheyenne_simulation/sim_data/LENS_piControl_HighRes_f09_g16_SOM_STATIC_20/ocn/hist/

ncrcat -O -v T_ML $d/LENS_piControl_HighRes_f09_g16_SOM_STATIC_20.ocn.h.monthly.{0081..0100}.nc all_rec.nc
for i in $(seq 1 12); do

    f2=clim_SST_$(printf "%02d" $i).nc
    
    ncra -O -F -d time,$i,,12 all_rec.nc $f2

done

ncrcat -O clim_SST_*.nc SST_SOM.nc

