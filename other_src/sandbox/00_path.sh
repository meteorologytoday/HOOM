#!/bin/bash

wkdir=`pwd`

gen_data_dir=$wkdir/gen_data
raw_data_dir=$wkdir/raw_data

#domain_file=$raw_data_dir/domain.lnd.fv4x5_gx3v7.091218.nc
domain_file=$raw_data_dir/fv45_remove_lnd.nc
#domain_file=$raw_data_dir/domain.ocn.gx3v7.120323.nc
#domain_file=$raw_data_dir/domain.ocn.gx1v6.090206.nc
record_file=$gen_data_dir/record_${1}.nc
#topo_file=$raw_data_dir/LENS_f45_g37_ocean_topog_gx1v6.nc
zdomain_file=$raw_data_dir/zdomain.nc
init_file=$gen_data_dir/init_aqp.nc
procs=${2}
