#!/bin/bash

for casename in lowres_STD_SOM lowres_SSM_SOM lowres_SSM_NK lowres_SSM_SOM_noQflux; do

    ncrcat -v ilev,TREFHT,VQ,VZ,VT,PRECC,PRECL,FSNT,FSNS,FLNT,FLNS,SHFLX,LHFLX $casename/atm/hist/$casename.cam.h0.00{01..20}-{01..12}.nc $casename.h0.nc
    #ncrcat -v ilev,T $casename/atm/hist/$casename.cam.h0.00{01..20}-{01..12}.nc $casename.h0.nc
    ncrcat -v ilev,TREFHT $casename/atm/hist/$casename.cam.h1.00{01..20}-01-01-00000.nc $casename.h1.nc

done
