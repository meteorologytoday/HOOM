#!/bin/bash

source 01_setting.sh

mkdir -p $out_dir
mkdir -p $tmp_dir

ptasks=6
set -x
# Doing average

cd $in_dir
for varname in TEMP SALT; do

    echo "Doing var: $varname"

    ((cnt=cnt%ptasks)); ((cnt++==0)) && wait

    
    year_rng=$( printf "{%04d..%04d}" $beg_y $end_y )
    eval "ncra -O -v ${varname} ${data_prefix_monthly}.${year_rng}-{01..12}.nc $out_dir/mean-${varname}.nc" & 

done

wait

cd $wkdir

echo "Done"
