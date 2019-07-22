#!/bin/bash

data_dir=./data
code_dir=./model-dev/other_src/stanfit
julia $code_dir/mk_docn_forcing_file.jl --input-file=$data_dir/forcing.gx3v7.nc --Qflux-varname=qflux_SOM --MLD-varname=h_SOM --domain-file=$data_dir/domain.ocn.gx3v7.120323.nc --output-file=$data_dir/docn_forcing.SOM_LENS.g37.nc
julia $code_dir/mk_docn_forcing_file.jl --input-file=$data_dir/forcing.gx3v7.nc --Qflux-varname=qflux_EntSOM15L --MLD-varname=h_EntSOM15L --domain-file=$data_dir/domain.ocn.gx3v7.120323.nc --output-file=$data_dir/docn_forcing.EntSOM15L_LENS.g37.nc
