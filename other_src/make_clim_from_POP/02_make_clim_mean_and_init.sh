#!/bin/bash

source 01_setting.sh

mkdir -p $out_dir
mkdir -p $tmp_dir

timeslab=$( printf "time,%d,%d,1" $(( ( beg_y - 1 ) * 12 + 1 + offset_month )) $(( ( end_y - 1 ) * 12 + 12 + offset_month )) )
init_time=$(( ( beg_y - 1 ) * 12 + 1 + offset_month ))

# Doing average and get initilization
for varname in HMXL TEMP SALT; do

    echo "Doing var: $varname"
    ncks -O -F -v $varname -d time,$init_time,$init_time $in_dir/${data_prefix_monthly}${varname}${data_suffix_monthly}.nc $out_dir/processed-init-${data_prefix_monthly}${varname}${data_suffix_monthly}.nc
    
    tmp_file=$tmp_dir/${varname}_yearly_avgs.nc
    rm $tmp_file

   ncra -O -F -d $timeslab $in_dir/${data_prefix_monthly}${varname}${data_suffix_monthly}.nc $out_dir/processed-avg-${data_prefix_monthly}${varname}${data_suffix_monthly}.nc &

done

wait

echo "Generate HMXL"
file=$out_dir/processed-${data_prefix_monthly}HMXL${data_suffix_monthly}.nc
ncks -O -F -v HMXL -d $timeslab $in_dir/${data_prefix_monthly}HMXL${data_suffix_monthly}.nc $file


echo "Generate SST"
# slice the surface
file=$out_dir/processed-${data_prefix_monthly}SST${data_suffix_monthly}.nc
ncks -O -F -v TEMP -d z_t,1,1 -d $timeslab $in_dir/${data_prefix_monthly}TEMP${data_suffix_monthly}.nc $file
ncwa -O -F -a z_t $file $file
ncrename -v TEMP,SST $file

echo "Generate SSS"
file=$out_dir/processed-${data_prefix_monthly}SSS${data_suffix_monthly}.nc
ncks -O -F -v SALT -d z_t,1,1 -d $timeslab $in_dir/${data_prefix_monthly}SALT${data_suffix_monthly}.nc $file
ncwa -O -F -a z_t $file $file
ncrename -v SALT,SSS $file


echo "Make HMXL into meters"
ncap2 -O -s "HMXL=HMXL/100.0" $out_dir/processed-avg-${data_prefix_monthly}HMXL${data_suffix_monthly}.nc $out_dir/processed-avg-${data_prefix_monthly}HMXL${data_suffix_monthly}.nc
ncap2 -O -s "HMXL=HMXL/100.0" $out_dir/processed-init-${data_prefix_monthly}HMXL${data_suffix_monthly}.nc $out_dir/processed-init-${data_prefix_monthly}HMXL${data_suffix_monthly}.nc
ncap2 -O -s "HMXL=HMXL/100.0" $out_dir/processed-${data_prefix_monthly}HMXL${data_suffix_monthly}.nc $out_dir/processed-${data_prefix_monthly}HMXL${data_suffix_monthly}.nc


