#!/bin/bash

for casename in lowres_STD_SOM lowres_SSM_SOM lowres_SSM_NK lowres_SSM_SOM_noQflux; do

    ncrcat -v ilev,VQ,VZ,VT,PRECC,PRECL $casename/atm/hist/$casename.cam.h0.00{11..20}-{01..12}.nc $casename.nc

done
