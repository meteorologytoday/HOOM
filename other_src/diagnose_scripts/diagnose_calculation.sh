#!/bin/bash

# Transform gx3v7 to fv45
if [ ! -f "$wgt_file" ]; then
    echo "Weight file \"$wgt_file\" does not exist, I am going to generate one..."
    julia -p 4  $coordtrans_scripts_path/generate_weight.jl --s-file=domain.ocn.gx3v7.120323.nc --d-file=domain.lnd.fv4x5_gx3v7.091218.nc --w-file=$wgt_file --s-mask-value=1.0 --d-mask-value=0.0
fi


# First, make a continuous file of atm/ocn output

echo "# Concat files..."
$diagnose_scripts_path/concat_files_atm.sh
$diagnose_scripts_path/concat_files_ocn.sh

# Second, transform grid
if [ ! -f "$ocn_trans_outputfile" ]; then
    julia $coordtrans_scripts_path/transform_data.jl --s-file=$ocn_outputfile --d-file=$ocn_trans_outputfile --w-file=$wgt_file --vars=T,MLD --x-dim=Nx --y-dim=Ny --t-dim=time 
fi


# Diagnose atm
echo "Diagnose atm..."
julia $analysis_path/atm_anomalies.jl --data-file=$atm_outputfile --domain-file=$atm_domain
julia $analysis_path/AO.jl --data-file-PSLA=$atm_outputfile_anomalies --domain-file=$atm_domain --EOF-file-AO=$AO_file


# Diagnose ocn in displaced pole grid
echo "Diagnose ocn..."
julia $analysis_path/SST_correlation.jl --data-file=$ocn_outputfile --domain-file=$ocn_domain --SST=T


# Transform processed data
julia $coordtrans_scripts_path/transform_data.jl --s-file=$ocn_outputfile --d-file=$ocn_trans_outputfile_SSTA --w-file=$wgt_file --vars=SSTA --x-dim=Nx --y-dim=Ny --z-dim=Nz --t-dim=time 
julia $coordtrans_scripts_path/transform_data.jl --s-file=$ocn_outputfile --d-file=$ocn_trans_outputfile_mstat --w-file=$wgt_file --vars=SSTAYYC,SSTAVAR --x-dim=Nx --y-dim=Ny --z-dim=Nz --t-dim=months


julia $analysis_path/PDO.jl --data-file-SSTA=$ocn_trans_outputfile_SSTA --domain-file=$atm_domain --EOF-file-PDO=$PDO_file


