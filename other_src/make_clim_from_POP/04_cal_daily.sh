#!/bin/bash

source 01_setting.sh

mkdir -p $tmp_dir
cd $tmp_dir

ptasks=6

varnames=( SST IFRAC )
new_varnames=( T_clim IFRAC_clim )

year_rng=$( printf "{%04d..%04d}" $beg_y $end_y )
set -x
for m in $(seq 1 12); do

    ((cnt=cnt%ptasks)); ((cnt++==0)) && wait
    {
        echo "# ncea month $m"

        for varname in "${varnames[@]}"; do
            eval "ncea -O -v $varname -h \
                $in_dir/${data_prefix_daily}.${year_rng}-$( printf '%02d' $m )-01.nc \
                tmp_${varname}_daily_$(printf '%02d' $m)-clim.nc "
        done
        
    } &
done

wait

for i in $( seq 0 $(( ${#varnames[@]} - 1 )) ); do

    varname="${varnames[$i]}"
    new_varname="${new_varnames[$i]}"
    fn=daily_clim_${new_varname}.nc

    ncrcat -3 -O -h tmp_${varname}_daily_{01..12}-clim.nc $fn
    ncrename -d nlat,nj -d nlon,ni -v ${varname},${new_varname} $fn
    
    mv $fn $out_dir
    
done


cd $wkdir
