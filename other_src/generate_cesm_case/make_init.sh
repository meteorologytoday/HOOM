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
    input-init-T-file
    input-init-S-file
    input-init-MLD-file
    input-topo-file
    output-clim-T-file
    output-clim-S-file
    output-init-T-file
    output-init-S-file
    output-init-MLD-file
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
if [ "$old_domain_file" == "$new_domain_file" ]; then
    wgt_file="X"
else
    wgt_file=$( basename $old_domain_file ".nc" )_$( basename $new_domain_file ".nc" ).nc
    
    wgt_dir="wgt_$( filename $ocn_domain )_to_$( filename $atm_domain )"
    $script_coordtrans_dir/generate_weight.sh    \
        --s-file=$old_domain_file                \
        --d-file=$new_domain_file                \
        --output-dir="$wgt_dir"                  \
        --s-mask-value=1.0                       \
        --d-mask-value=1.0

fi


# Convert 3D variable: TEMP, SALT
data_files=(
    TEMP $input_clim_T_file $output_clim_T_file
    SALT $input_clim_S_file $output_clim_S_file
    TEMP $input_init_T_file $output_init_T_file
    SALT $input_init_S_file $output_init_S_file
)

for i in $( seq 1 $(( ${#data_files[@]} / 3))); do
    varname=${data_files[$((3*(i-1)))]}
    data_file=${data_files[$((3*(i-1)+1))]}
    new_data_file=${data_files[$((3*(i-1)+2))]}

    echo "[$varname] $data_file => $new_data_file"

    tmp1=$tmp_dir/${label}_$( basename $data_file ".nc" ).clim.nc
    tmp2=$tmp_dir/${label}_${resolution}_$( basename $data_file ".nc" ).new-domain.nc

    if [ ! -f $new_data_file ]; then
        
        echo "Transforming variable: $varname"
        
        ncwa -O -a time $data_file $tmp1
        ncks -O -3 $tmp1 $tmp1
        ncrename -d nlat,Ny -d nlon,Nx -d z_t,Nz $tmp1
        ncks -O -4 $tmp1 $tmp1

        # Horizontal resolution
        if [ "$wgt_file" != "X" ]; then
            julia $script_coordtrans_dir/transform_data.jl --s-file=$tmp1 --d-file=$tmp2 --w-file=${wgt_dir}/wgt.bilinear.nc --vars=$varname --x-dim=Nx --y-dim=Ny --z-dim=Nz --algo=ESMF
        else
            mv $tmp1 $tmp2
        fi

        julia $script_coordtrans_dir/convert_z.jl $tmp2 $new_data_file $varname
        
#        rm -f $tmp1 $tmp2
    fi
done


# Convert 2D variable: MLD, TOPO
data_files=(
    HMXL   $input_init_MLD_file $output_init_MLD_file
    depth  $input_topo_file $output_topo_file
)

for i in $( seq 1 $(( ${#data_files[@]} / 3))); do
    varname=${data_files[$((3*(i-1)))]}
    data_file=${data_files[$((3*(i-1)+1))]}
    new_data_file=${data_files[$((3*(i-1)+2))]}

    tmp=$tmp_dir/${label}_$( basename ${data_file} ".nc" ).tmp.nc
    if [ ! -f "$new_data_file" ]; then

        echo "Transforming variable: $varname"

        ncks -O -3 $data_file $tmp
        ncrename -d ni,Nx -d nj,Ny $tmp 
        ncks -O -4 $tmp $tmp
     

        if [ "$wgt_file" != "X" ]; then
            julia $script_coordtrans_dir/transform_data.jl --s-file=$tmp --d-file=$new_data_file --w-file=${wgt_dir}/wgt.bilinear.nc --vars=$varname --x-dim=Nx --y-dim=Ny --algo=ESMF
        else
            mv $tmp $new_data_file
        fi
    fi
done

# Make z-coordinate file
if [ ! -f $output_zdomain_file ]; then
    julia $script_coordtrans_dir/SSM_mk_zdomain.jl $output_zdomain_file
fi



