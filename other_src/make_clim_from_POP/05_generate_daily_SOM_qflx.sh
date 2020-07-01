#!/bin/bash

source 01_setting.sh

julia $wkdir/mk_qflux_forcing/interpolate_monthly_to_daily.jl \
    --domain-file=${domain_file}                       \
    --input-file=init_guess_qflux_T.nc                 \
    --varname=qflux_T            \
    --time=mid                   \
    --output-file=$out_dir/SOM_qflx_T.nc
