#!/bin/bash

set -x
source 00_path.sh

analysis_file=$gen_data_dir/analysis_${1}.nc
analysis_file2=$gen_data_dir/analysis_weighted_${1}.nc

cp $record_file $analysis_file

ncks -O -3 $domain_file tmp.nc
ncrename -d ni,Nx -d nj,Ny tmp.nc

ncks -A -v area tmp.nc $analysis_file
ncks -A -v mask $init_file $analysis_file


rm tmp.nc

ncwa -O -b -a Nx,Ny -B 'mask == 1' -w area -v TEMP,dTEMPdt,TFLUX_bot,SALT,dSALTdt,SFLUX_bot,SFLUX_top,TSAS_clim,SSAS_clim,TFLUX_DIV_implied,SFLUX_DIV_implied $analysis_file $analysis_file2
#ncap2 -O -s 'TFLUX_DIV_implied' tmp.nc $analysis_file2

#rm tmp.nc
