#!/bin/bash

export script_root_dir=$(dirname $0)
export script_coordtrans_dir=$script_root_dir/../../CoordTrans
export s_domain=domain.ocn.gx1v6.090206.nc
export d_domain=domain.ocn.gx3v7.120323.nc
export wgt_file=wgt_gx1v6_to_gx3v7.nc

i_fmt="b.e11.B1850C5CN.f09_g16.005.pop.h.%s.100001-109912.nc"
o_fmt="transformed_b.e11.B1850C5CN.f45_g37.005.pop.h.%s.100001-109912.nc"

vars_time=( HMXL SHF SST )
vars_avg=( TEMP )

if [ -f "$wgt_file" ]; then
    echo "Weight file \"$wgt_file\" already exists. No need to create a new one."
else
    echo "Weight file \"$wgt_file\" does not exist, I am going to generate one..."
    julia -p 4  $script_coordtrans_dir/generate_weight.jl --s-file=$s_domain --d-file=$d_domain --w-file=$wgt_file --s-mask-value=1.0 --d-mask-value=1.0
fi

# transform
for var in "${vars_time[@]}"; do
    i_file=$( printf $i_fmt  $var )
    o_file=$( printf $o_fmt  $var )

    if [ -f "$o_file" ]; then
        echo "$o_file already exist. No need to do anything."
        
    else
        echo "$o_file does not exist. Now we produce it..."
        julia $script_coordtrans_dir/transform_data.jl --s-file=$i_file --d-file=$o_file --w-file=$wgt_file --vars=$var --x-dim=nlon --y-dim=nlat --z-dim=z_t --t-dim=time
    fi
done

 # transform
for var in "${vars_avg[@]}"; do
    i_file=$( printf $i_fmt  $var )
    o_file=$( printf $o_fmt  $var )

    tmp=$i_file.tmp

    if [ -f "$o_file" ]; then
        echo "$o_file already exist. No need to do anything."
        
    else
        echo "$o_file does not exist. Now we produce it..."
        ncra -O $i_file $tmp
        julia $script_coordtrans_dir/transform_data.jl --s-file=$tmp --d-file=$o_file --w-file=$wgt_file --vars=$var --x-dim=nlon --y-dim=nlat --z-dim=z_t --t-dim=time
        rm $tmp

        ncks -A -v z_w_top,z_w_bot $i_file $o_file
    fi
done 
