#!/bin/bash

lopts=(res casenames sim-data-dir diag-data-dir graph-dir atm-domain ocn-domain colors)

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

if [ -f flag_plot_mc ]; then
    echo "foldeer: $script_plot_dir"
    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=ocn_concat_rg.nc --varname=PDO --colors="$colors"
    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=atm_analysis5_AO.nc --varname=AO --colors="$colors"
    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=ocn_concat_rg.nc --varname=EN34 --colors="$colors"

    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=atm_analysis2.nc --varname=TREFHT_GLB --ylabel="Temperature [K]" --mavg=1 --extra-title="Yearly average" --colors="$colors"
    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=atm_analysis2.nc --varname=TREFHT_LND --ylabel="Temperature [K]" --mavg=1 --extra-title="Yearly average" --colors="$colors"
    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=atm_analysis2.nc --varname=TREFHT_OCN --ylabel="Temperature [K]" --mavg=1 --extra-title="Yearly average" --colors="$colors"
    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=ice_trends.nc --varname=total_ice_volume --ylabel="Ice volume [ \$ \\mathrm{km^3} \$ ]" --yscale="1e9" --colors="$colors"
    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=ice_trends.nc --varname=total_ice_area --ylabel="Ice area fraction [ % ]" --yscale="1e-2" --colors="$colors"

    python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=atm_analysis3_energy.nc --varname=IAET_mean --ylabel="Energy flux [ \$ \\times 10^{15} \\, \\mathrm{W} \$ ]" --yscale="1e15" --domain-file=$atm_domain --colors="$colors"
    python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=atm_analysis4_precip.nc --varname=total_precip_ZM --ylabel="Precipitation [ \$ \\mathrm{mm} \, \\mathrm{yr}^{-1} \$ ]" --yscale="3.171e-11" --domain-file=$atm_domain --colors="$colors"
    python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casenames --data-file=atm_analysis4_precip.nc --varname=total_precip_ZVAR --ylabel="Precipitation [ \$ \\mathrm{mm} \, \\mathrm{yr}^{-1} \$ ]" --yscale="3.171e-11" --domain-file=$atm_domain --colors="$colors"

fi
