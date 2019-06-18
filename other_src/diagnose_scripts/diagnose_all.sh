#!/bin/bash

export casename=$1
export archive_path=$2
export nc_output_dir=$3
export diagnose_output_dir=$4
export beg_year=$(printf "%04d" $5)
export end_year=$(printf "%04d" $6)

export atm_hist_path=$archive_path/atm/hist
export ocn_hist_path=$archive_path/ocn

export ocn_outputfile=$nc_output_dir/$casename.ocn.h.monthly.nc
export year_stamp=$(printf "%s-%s" $beg_year $end_year)

export atm_domain=domain.lnd.fv4x5_gx3v7.091218.nc
export ocn_domain=domain.ocn.gx3v7.120323.nc

export wpath=`pwd`
export diagnose_scripts_path=$(dirname $0)
export analysis_path=$diagnose_scripts_path/analysis
export coordtrans_scripts_path=$diagnose_scripts_path/../CoordTrans
export wgt_file=wgt_gx3v7_to_fv4x5.nc
export PDO_file=PDO_EOFs_fv45.nc
export AO_file=AO_EOFs_fv45.nc


for dir_path in $nc_output_dir $diagnose_output_dir ; do

    echo "Checking path: $dir_path"
    if [ ! -d $dir_path ]; then
        mkdir -p $dir_path
    fi

done


# filenames
export atm_outputfile=$nc_output_dir/$casename.atm.h0.${year_stamp}.nc
export atm_outputfile_anomalies=$nc_output_dir/atm_anomalies.nc

export ocn_outputfile=$nc_output_dir/$casename.ocn.h.monthly.${year_stamp}.nc
export ocn_trans_outputfile=$nc_output_dir/$casename.ocn.h.monthly.transformed.${year_stamp}.nc

export ocn_trans_outputfile_SSTA=$nc_output_dir/$casename.ocn.h.monthly.SSTA.transformed.${year_stamp}.nc
export ocn_trans_outputfile_mstat=$nc_output_dir/$casename.ocn.h.monthly.mstat.transformed.${year_stamp}.nc

if [ ! -f flag_nocal ]; then
    $diagnose_scripts_path/diagnose_calculation.sh
fi

if [ ! -f flag_noplot ]; then
    $diagnose_scripts_path/diagnose_plot.sh
fi
