#!/bin/bash

wk_dir=$( dirname $0 )
script_coordtrans_dir=$wk_dir/../CoordTrans
tmp_dir=tmp

echo "wk_dir: $wk_dir"

lopts=(output-dir label data-clim-T-file data-clim-S-file old-domain-file new-domain-file topo-file T-unit)

options=$(getopt -o '' --long $(printf "%s:," "${lopts[@]}") -- "$@")
[ $? -eq 0 ] || { 
    echo "Incorrect options provided"
    exit 1
}
eval set -- "$options"


while true; do
    for lopt in "${lopts[@]}"; do
        eval "if [ \"\$1\" == \"--$lopt\" ]; then shift; export ${lopt//-/_}=\"\$1\"; shift; break; fi"
    done

    if [ "$1" == -- ]; then
        shift;
        break;
    fi
done

echo "Received parameters: "
for lopt in "${lopts[@]}"; do
    llopt=${lopt//-/_}
    eval "echo \"- $llopt=\$$llopt\""
done

IFS=', ' read -r -a data_clim_T_file <<< "$data_clim_T_file"
IFS=', ' read -r -a data_clim_S_file <<< "$data_clim_S_file"
IFS=', ' read -r -a topo_file        <<< "$topo_file"


case_settings=(
    SOM                      make_SOM_init.jl
    ESOM                     make_ESOM_init.jl
    MLMML_restricted_climate make_MLMML_init_restricted_climate.jl
    MLMML_unrestricted       make_MLMML_init_unrestricted.jl
)

case_settings=(
    ESOM                     make_ESOM_init.jl
)


mkdir -p $output_dir
mkdir -p $tmp_dir

# First generate correct transformed  coordinate files

wgt_file=$( basename $old_domain_file ".nc" )_$( basename $new_domain_file ".nc" ).nc
new_data_clim_T_file=${label}_$( basename $data_clim_T_file ".nc" )

if [ ! -f "$wgt_file" ]; then
    echo "Weight file \"$wgt_file\" does not exist, I am going to generate one..."
    julia -p 4  $script_coordtrans_dir/generate_weight.jl --s-file=$old_domain_file --d-file=$new_domain_file --w-file=$wgt_file --s-mask-value=1.0 --d-mask-value=1.0
fi

data_files=(
    ${data_clim_T_file[@]}
    ${data_clim_S_file[@]}
)

new_data_files=()
for i in $( seq 1 $(( ${#data_files[@]} / 2))); do
    data_file=${data_files[$((2*(i-1)))]}
    varname=${data_files[$((2*(i-1)+1))]}

    tmp1=$tmp_dir/${label}_$( basename $data_file ".nc" ).clim.nc
    tmp2=$tmp_dir/${label}_$( basename $data_file ".nc" ).new-domain.nc
    new_data_file=$output_dir/${label}_$( basename $data_file ".nc" ).nc
    new_data_files+=( "$new_data_file,$varname" )

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
tmp=$tmp_dir/${label}_$( basename ${topo_file[0]} ".nc" ).tmp.nc
new_topo_file=$output_dir/${label}_$( basename ${topo_file[0]} ".nc" ).nc
if [ ! -f $new_topo_file ]; then

    ncks -O -3 ${topo_file[0]} $tmp
    ncrename -d ni,Nx -d nj,Ny $tmp 
    ncks -O -4 $tmp $tmp
    
    julia $script_coordtrans_dir/transform_data.jl --s-file=$tmp --d-file=$new_topo_file --w-file=$wgt_file --vars=${topo_file[1]} --x-dim=Nx --y-dim=Ny 

    rm -f $tmp
fi

# Make z-coordinate file
tmp_zdomain=$tmp_dir/z_domain.nc
julia $script_coordtrans_dir/SSM_mk_zdomain.jl $tmp_zdomain


# Make init files
for i in $(seq 1 $((${#case_settings[@]}/2))); do
    casename=${case_settings[$((2*(i-1)))]}
    gen_code=${case_settings[$((2*(i-1)+1))]}
    printf "[%s] => [%s]\n" $casename $gen_code

    output_file=$output_dir/init_${label}_${casename}.nc

    julia $wk_dir/init_code/$gen_code --output-file=$output_file --data-clim-T-file=${new_data_files[0]} --data-clim-S-file=${new_data_files[1]} --domain-file=$new_domain_file,ni,nj,mask --domain-z-file=$tmp_zdomain,zs --topo-file=$new_topo_file,depth --T-unit=$T_unit

done

