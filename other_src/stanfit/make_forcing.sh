#!/bin/bash

data_dir=../../data

julia mk_docn_forcing_file.jl --Qflux-file=$data_dir/SOM_fixed_MLD_LENS.g37_c2_s1000_w100.nc --Qflux-varname=Q_mean --domain-file=$data_dir/domain.ocn.gx3v7.120323.nc --output-file=$data_dir/docn_forcing.SOM_fixed_MLD_LENS.g37.nc
