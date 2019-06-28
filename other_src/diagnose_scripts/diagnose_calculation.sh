#!/bin/bash


if [ ! -f flag_nodiag_atm ] ; then

    # Diagnose atm
    echo "Diagnose atm..."
    julia $script_analysis_dir/atm_anomalies.jl --data-file=$atm_concat --domain-file=$atm_domain --output-file=$atm_analysis1
    julia $script_analysis_dir/atm_temperature.jl --data-file=$atm_concat --domain-file=$atm_domain --output-file=$atm_analysis2
    julia $script_analysis_dir/AO.jl --data-file-PSLA=$atm_analysis1 --domain-file=$atm_domain --EOF-file-AO=$AO_file
fi

if [ ! -f flag_nodiag_ocn ] ; then


    # Transform processed data
#    julia $coordtrans_scripts_path/transform_data.jl --s-file=$ocn_concat --d-file=$ocn_concat_rg_anomalies --w-file=$wgt_file --vars=SSTA --x-dim=Nx --y-dim=Ny --z-dim=Nz --t-dim=time 
 #   julia $coordtrans_scripts_path/transform_data.jl --s-file=$ocn_concat --d-file=$ocn_concat_rg_mstat --w-file=$wgt_file --vars=SSTAYYC,SSTAVAR --x-dim=Nx --y-dim=Ny --z-dim=Nz --t-dim=months

    # Diagnose ocn 
    echo "Diagnose ocn..."
    julia $script_analysis_dir/SST_correlation.jl --data-file=$ocn_concat_rg --domain-file=$atm_domain --SST=T
    julia $script_analysis_dir/PDO.jl --data-file-SSTA=$ocn_concat_rg --domain-file=$atm_domain --EOF-file-PDO=$PDO_file
    julia $script_analysis_dir/EN34.jl --data-file-SSTA=$ocn_concat_rg --domain-file=$atm_domain 

fi

if [ ! -f flag_nodiag_ice ] ; then

    echo "Diagnose ice..."
    julia $script_analysis_dir/ice.jl --data-file=$ice_concat_rg --domain-file=$atm_domain --output-file=$ice_analysis1

fi
