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

year_stamp=$(printf "%s-%s" $beg_year $end_year)


for dir_path in $nc_output_dir $diagnose_output_dir ; do

    echo "Checking path: $dir_path"
    if [ ! -d $dir_path ]; then
        mkdir -p $dir_path
    fi

done

wpath=`pwd`
diagnose_scripts_path=$(dirname $0)
coordtrans_scripts_path=$diagnose_scripts_path/../CoordTrans
wgt_file=wgt_gx3v7_to_fv4x5.nc

atm_domain=domain.lnd.fv4x5_gx3v7.091218.nc
ocn_domain=domain.ocn.gx3v7.120323.nc

# Transform gx3v7 to fv45
if [ ! -f "$wgt_file" ]; then
    echo "Weight file \"$wgt_file\" does not exist, I am going to generate one..."
    julia -p 4  $coordtrans_scripts_path/generate_weight.jl --s-file=domain.ocn.gx3v7.120323.nc --d-file=domain.lnd.fv4x5_gx3v7.091218.nc --w-file=$wgt_file --s-mask-value=1.0 --d-mask-value=0.0
fi



# First, make a continuous file of atm/ocn output


export atm_outputfile=$nc_output_dir/$casename.atm.h0.${year_stamp}.nc
export ocn_outputfile=$nc_output_dir/$casename.ocn.h.monthly.${year_stamp}.nc

export ocn_trans_outputfile=$nc_output_dir/$casename.ocn.h.monthly.transformed.${year_stamp}.nc

echo "# Concat files..."
$diagnose_scripts_path/concat_files_atm.sh
$diagnose_scripts_path/concat_files_ocn.sh


# Second, transform grid
if [ ! -f "$ocn_trans_outputfile" ]; then
    julia $coordtrans_scripts_path/transform_data.jl --s-file=$ocn_outputfile --d-file=$ocn_trans_outputfile --w-file=$wgt_file --vars=T,MLD --x-dim=Nx --y-dim=Ny --t-dim=time 
fi

# Diagnose atm
echo "Diagnose atm..."
python3 $diagnose_scripts_path/plot_SST.py --data-file=$atm_outputfile --domain-file=$atm_domain --output-dir=$diagnose_output_dir

# Diagnose ocn
echo "Diagnose ocn..."
#julia $diagnose_scripts_path/SST_correlation.jl --data-file=$ocn_outputfile --domain-file=$ocn_domain --SST=T



# Transform processed data

export ocn_trans_outputfile_SSTA=$nc_output_dir/$casename.ocn.h.monthly.SSTA.transformed.${year_stamp}.nc
export ocn_trans_outputfile_mstat=$nc_output_dir/$casename.ocn.h.monthly.mstat.transformed.${year_stamp}.nc

#julia $coordtrans_scripts_path/transform_data.jl --s-file=$ocn_outputfile --d-file=$ocn_trans_outputfile_SSTA --w-file=$wgt_file --vars=SSTA --x-dim=Nx --y-dim=Ny --t-dim=time 

#julia $coordtrans_scripts_path/transform_data.jl --s-file=$ocn_outputfile --d-file=$ocn_trans_outputfile_mstat --w-file=$wgt_file --vars=SSTAYYC,SSTAVAR --x-dim=Nx --y-dim=Ny --t-dim=months

 # Plotting 
python3 $diagnose_scripts_path/plot_ocean_diagnose.py --data-file-SSTAYYC=$ocn_trans_outputfile_mstat --data-file-SSTAVAR=$ocn_trans_outputfile_mstat --domain-file=$atm_domain --output-dir=$diagnose_output_dir


