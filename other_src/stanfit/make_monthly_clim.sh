#!/bin/bash

mkdir -p tmp


for i in $( seq 1 12 ); do

    ncra -O -d time,$i,,12 -v HMXL transformed_b.e11.B1850C5CN.f45_g37.005.pop.h.HMXL.100001-109912.nc tmp/$( printf "%02d" $i ).nc

done

ncrcat -O -n 12,2,1 tmp/01.nc HMXL_clim.nc
ncap2 -v -O -s "HMXL=HMXL/100" HMXL_clim.nc HMXL_clim_meters.nc

rm -rf tmp



