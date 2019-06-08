#!/bin/bash

export casename=$1
export archive_path=$2
export nc_output_dir=$3
export diagnose_output_dir=$4
export beg_year=$(printf "%04d" $5)
export end_year=$(printf "%04d" $6)

export atm_hist_path=$archive_path/atm/hist
export ocn_hist_path=$archive_path/ocn/hist

for dir_path in $nc_output_dir $diagnose_output_dir ; do

    echo "Checking path: $dir_path"
    if [ ! -d $dir_path ]; then
        mkdir -p $dir_path
    fi

done

wpath=`pwd`
diagnose_scripts_path=$(dirname $0)



# First, make a continuous file of atm/ocn output
$diagnose_scripts_path/concat_files.sh


# Diagnose atm
python3 $diagnose_scripts_path/plot_SST.py $nc_output_dir/$casename.h0.nc $nc_output_dir/$casename.h0.nc


# Diagnose ocn



