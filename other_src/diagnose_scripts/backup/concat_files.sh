#!/bin/bash


out_dir=extract_nc_files

if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
fi


for casename in lowres_STD_SOM lowres_SSM_SOM lowres_SSM_NK lowres_SSM_SOM_noQflux; do
    in_dir="../$casename/atm/hist"

    printf "Processing %s...\n" $casename

    ncrcat -O -v ilev,PSL,V,TREFHT,VQ,VZ,VT,PRECC,PRECL,FSNT,FSNS,FLNT,FLNS,SHFLX,LHFLX $in_dir/$casename.cam.h0.00{01..20}-{01..12}.nc $out_dir/$casename.h0.nc
    #ncrcat -v ilev,T $casename/atm/hist/$casename.cam.h0.00{01..20}-{01..12}.nc $casename.h0.nc
    ncrcat -O -v ilev,TREFHT $in_dir/$casename.cam.h1.00{01..20}-01-01-00000.nc $out_dir/$casename.h1.nc

done
