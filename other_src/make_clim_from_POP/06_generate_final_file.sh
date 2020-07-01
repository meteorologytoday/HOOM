#!/bin/bash

source 01_setting.sh

stamp=$( date +'%Y%m%d' )

ConH=$out_dir/docn_forcing.${stamp}.g16.daily.init.ConH.nc
rm $ConH

rm $daily_clim_file
echo "Creating $daily_clim_file"

merged_file=$tmp_dir/merged.nc

rm $merged_file
for varname in MLD S_clim T_clim IFRAC_clim; do
    echo "Appending $varname ..."
    ncks -A -v $varname $out_dir/daily_clim_${varname}.nc $merged_file
done

julia mk_qflux_forcing/mk_docn_forcing_file_daily.jl \
    --input-Qflx-file=""                     \
    --input-clim-file=$merged_file       \
    --input-MLD-file=$merged_file        \
    --Qflx_T-varname=qflux_T                 \
    --domain-file=$domain_file               \
    --output-file=$ConH                      \
    --MLD-average
