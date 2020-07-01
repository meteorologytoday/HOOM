#!/bin/bash

source 01_setting.sh

mkdir -p $tmp_dir
cd $tmp_dir

ptasks=6

year_rng=$( printf "[%04d-%04d]" $beg_year $end_year )
for m in $(seq 1 12); do

    ((cnt=cnt%ptasks)); ((cnt++==0)) && wait
    {
        echo "# ncea month $m"

        for varname in SST IFRAC; do
            ncea -O -v $varname -h \
                "$in_dir/${data_prefix_daily}.${year_rng}-$( printf '%02d' $m ).nc" \
                tmp_${varname}_daily_$(printf "%02d" $m)-clim.nc
        done
        
    } &
done

wait

set -x
for varname in SST IFRAC; do
    ncrcat -3 -O -h tmp_${varname}_daily_[01-12]-clim.nc tmp_${varname}_daily_clim.nc
done

mv tmp_SST_daily_clim.nc  daily_clim_T_clim.nc

ncrename -d nlat,nj -d nlon,ni daily_clim_T_clim.nc
ncrename -v SST,T_clim -v IFRAC,IFRAC_clim daily_clim_T_clim.nc

mv daily_clim_T_clim.nc $out_dir

cd $wkdir
