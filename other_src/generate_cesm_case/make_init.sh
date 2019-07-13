#!/bin/bash

#wk_dir=$( dirname $0 )
script_coordtrans_dir=$wk_dir/../CoordTrans
tmp_dir=tmp

echo "wk_dir: $wk_dir"

lopts=(
    output-dir
    label
    input-clim-T-file
    input-clim-S-file
    input-topo-file
    output-clim-T-file
    output-clim-S-file
    output-topo-file
    old-domain-file
    new-domain-file
    T-unit
    output-zdomain-file
)

source $wk_dir/getopt_helper.sh

mkdir -p $output_dir
mkdir -p $tmp_dir

# First generate correct transformed  coordinate files

wgt_file=$( basename $old_domain_file ".nc" )_$( basename $new_domain_file ".nc" ).nc

if [ ! -f "$wgt_file" ]; then
    echo "Weight file \"$wgt_file\" does not exist, I am going to generate one..."
    julia -p 4  $script_coordtrans_dir/generate_weight.jl --s-file=$old_domain_file --d-file=$new_domain_file --w-file=$wgt_file --s-mask-value=1.0 --d-mask-value=1.0
fi

data_files=(
    TEMP $input_clim_T_file $output_clim_T_file
    SALT $input_clim_S_file $output_clim_S_file
)

for i in $( seq 1 $(( ${#data_files[@]} / 3))); do
    varname=${data_files[$((3*(i-1)))]}
    data_file=${data_files[$((3*(i-1)+1))]}
    new_data_file=${data_files[$((3*(i-1)+2))]}

    echo "[$varname] $data_file => $new_data_file"

    tmp1=$tmp_dir/${label}_$( basename $data_file ".nc" ).clim.nc
    tmp2=$tmp_dir/${label}_$( basename $data_file ".nc" ).new-domain.nc

    if [ ! -f $new_data_file ]; then
        
        ncra -O -d time,1 $data_file $tmp1
        ncks -O -3 $tmp1 $tmp1
        ncrename -d nlat,Ny -d nlon,Nx -d z_t,Nz $tmp1
        ncks -O -4 $tmp1 $tmp1

        # Horizontal resolution
        julia $script_coordtrans_dir/transform_data.jl --s-file=$tmp1 --d-file=$tmp2 --w-file=$wgt_file --vars=$varname --x-dim=Nx --y-dim=Ny --z-dim=Nz --t-dim=time 
        julia $script_coordtrans_dir/convert_z.jl $tmp2 $new_data_file $varname
        
        rm -f $tmp1 $tmp2
    fi
done


# Make topography
tmp=$tmp_dir/${label}_$( basename ${input_topo_file} ".nc" ).tmp.nc
if [ ! -f $output_topo_file ]; then

    ncks -O -3 ${topo_file[0]} $tmp
    ncrename -d ni,Nx -d nj,Ny $tmp 
    ncks -O -4 $tmp $tmp
    
    julia $script_coordtrans_dir/transform_data.jl --s-file=$tmp --d-file=$output_topo_file --w-file=$wgt_file --vars=depth --x-dim=Nx --y-dim=Ny 

fi

# Make z-coordinate file
if [ ! -f $output_zdomain_file ]; then
    julia $script_coordtrans_dir/SSM_mk_zdomain.jl $output_zdomain_file
fi



