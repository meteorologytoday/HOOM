#!/bin/bash

source 01_setting.sh

mkdir -p $tmp_dir
mkdir -p $out_dir
cd $tmp_dir

ptasks=6



year_rng=$( printf "{%04d..%04d}" $beg_y $end_y )

set -x
for m in $(seq 1 12); do

    ((cnt=cnt%ptasks)); ((cnt++==0)) && wait
    {
        echo "# ncea month $m"

        for varname in SALT HMXL; do
            eval "ncea -O -v $varname -h \
                $in_dir/${data_prefix_monthly}.${year_rng}-$( printf '%02d' $m ).nc \
                tmp_${varname}_monthly_$(printf '%02d' $m)-clim.nc"
        done
    } &
done

wait

echo "concat monthly files"
set -x
for varname in SALT HMXL; do
    ncrcat -3 -O -h -v $varname tmp_${varname}_monthly_{01..12}-clim.nc tmp_${varname}_monthly_clim.nc
done

monthly_clim_file=$tmp_dir/monthly_clim.nc

# create clim file using SSS    
#ncks -O -3 -F -h -v SALT -d z_t,1,1 tmp_SALT_monthly_{01..12}-clim.nc tmp_SSS_monthly_clim.nc
ncwa -a z_t tmp_SALT_monthly_clim.nc $monthly_clim_file
ncrename -v SALT,SSS $monthly_clim_file

# Deal with mixed-layer
ncap2 -O -s 'HMXL=HMXL/100.0;' tmp_HMXL_monthly_clim.nc tmp_HMXL_monthly_clim.nc
ncks -A -v HMXL  tmp_HMXL_monthly_clim.nc $monthly_clim_file


ncrename -d nlat,nj -d nlon,ni $monthly_clim_file
ncrename -v SSS,S_clim -v HMXL,MLD $monthly_clim_file

for varname in S_clim MLD; do
    julia $wkdir/mk_qflux_forcing/interpolate_monthly_to_daily.jl \
        --domain-file=${domain_file} \
        --input-file=${monthly_clim_file} \
        --varname=${varname}         \
        --time=mid                   \
        --output-file=daily_clim_${varname}.nc
done


mv monthly_clim.nc $out_dir
mv daily_clim_*.nc $out_dir

cd $wkdir
