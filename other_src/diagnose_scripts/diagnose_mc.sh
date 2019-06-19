#!/bin/bash

lopts=(res casenames sim-data-dir diag-data-dir graph-dir atm-domain ocn-domain)

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


export wdir=`pwd`
export script_root_dir=$(dirname $0)
script_plot_dir=$script_root_dir/plot

mkdir -p $graph_dir

if [ ! -f flag_noplot ]; then
    echo "foldeer: $script_plot_dir"
    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=ocn_concat_rg.nc --varname=PDO
    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=ocn_concat_rg.nc --varname=EN34

    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=atm_analysis2.nc --varname=TREFHT_GLB --normalize=no --ylabel="Temperature [K]"
    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=atm_analysis2.nc --varname=TREFHT_LND --normalize=no --ylabel="Temperature [K]"
    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=atm_analysis2.nc --varname=TREFHT_OCN --normalize=no --ylabel="Temperature [K]"
fi
