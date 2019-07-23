#!/bin/bash

export wk_dir=$( dirname $0 )
script_coordtrans_dir=$wk_dir/../CoordTrans

lopts=(
    input-file
    output-dir
)

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

old_domain_file=/seley/tienyiah/CESM_domains/domain.ocn.gx3v7.120323.nc
new_domain_file=/seley/tienyiah/CESM_domains/domain.lnd.fv4x5_gx3v7.091218.nc


# First generate correct transformed  coordinate files
wgt_file=$( basename $old_domain_file ".nc" )_$( basename $new_domain_file ".nc" ).nc

if [ ! -f "$wgt_file" ]; then
    echo "Weight file \"$wgt_file\" does not exist, I am going to generate one..."
    julia -p 4 $script_coordtrans_dir/generate_weight.jl --s-file=$old_domain_file --d-file=$new_domain_file --w-file=$wgt_file --s-mask-value=1.0 --d-mask-value=0.0
fi

new_input_file="trans_$( basename $input_file )"

if [ ! -f "$new_input_file" ]; then 
    julia $script_coordtrans_dir/transform_data.jl --w-file=$wgt_file --s-file=$input_file --d-file=$new_input_file --x-dim=Nx --y-dim=Ny --z-dim=Nz --t-dim=time
fi

for i in $( seq 1 12 ); do

    python3 $wk_dir/../Qflux/plot_scripts/plot_contourf.py \
        --data-file=$new_input_file    \
        --domain-file=$new_domain_file \
        --output-dir=$output_dir       \
        --casename=EntOM-xM            \
        --varname=qflux_EntSOM15L      \
        --colormap=cmocean_balance     \
        --cmin="-200"                  \
        --cmax="200"                   \
        --clevs=20                     \
        --idx-t=$(( $i - 1 ))          \
        --clabel="[ W / m^2 ]"         \
        --title="[EntOM-xM] Qflux of month $(printf "%02d" $i)"\
        --extra-filename="$( printf "_%02d" $i )"

    python3 $wk_dir/../Qflux/plot_scripts/plot_contourf.py \
        --data-file=$new_input_file    \
        --domain-file=$new_domain_file \
        --output-dir=$output_dir       \
        --casename=EntOM-xM            \
        --varname=eflux_EntSOM15L      \
        --colormap=cmocean_balance     \
        --cmin="-200"                  \
        --cmax="200"                   \
        --clevs=20                     \
        --idx-t=$(( $i - 1 ))          \
        --clabel="[ W / m^2 ]"         \
        --title="[EntOM-xM] Eflux of month $(printf "%02d" $i)"\
        --extra-filename="$( printf "_%02d" $i )"

    python3 $wk_dir/../Qflux/plot_scripts/plot_contourf.py \
        --data-file=$new_input_file    \
        --domain-file=$new_domain_file \
        --output-dir=$output_dir       \
        --casename=EntOM-xM            \
        --varname=h_EntSOM15L          \
        --colormap=cmocean_deep        \
        --cmin="0"                     \
        --cmax="200"                   \
        --clevs=20                     \
        --idx-t=$(( $i - 1 ))          \
        --clabel="[ m ]"               \
        --title="[EntOM-xM] MLD of month $(printf "%02d" $i)"\
        --extra-filename="$( printf "_%02d" $i )"

    python3 $wk_dir/../Qflux/plot_scripts/plot_contourf.py \
        --data-file=$new_input_file    \
        --domain-file=$new_domain_file \
        --output-dir=$output_dir       \
        --casename=SOM                 \
        --varname=qflux_SOM            \
        --colormap=cmocean_balance     \
        --cmin="-200"                  \
        --cmax="200"                   \
        --clevs=20                     \
        --idx-t=$(( $i - 1 ))          \
        --clabel="[ W / m^2 ]"         \
        --title="[SOM] Qflux of month $(printf "%02d" $i)"\
        --extra-filename="$( printf "_%02d" $i )"

done


python3 $wk_dir/../Qflux/plot_scripts/plot_contourf.py \
    --data-file=$new_input_file    \
    --domain-file=$new_domain_file \
    --output-dir=$output_dir       \
    --casename=SOM                 \
    --varname=h_SOM                \
    --colormap=cmocean_deep        \
    --cmin="0"                     \
    --cmax="200"                   \
    --clevs=20                     \
    --idx-t=0                      \
    --clabel="[ m ]"               \
    --title="[SOM] MLD"            

