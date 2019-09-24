#!/bin/bash

lopts=(
    res
    label
    casenames
    legends
    sim-data-dir
    diag-data-dir
    graph-data-dir
    atm-domain
    ocn-domain
    colors
    linestyles
    diag-beg-year
    diag-end-year
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


export wdir=`pwd`
export script_root_dir=$(dirname $0)
script_plot_dir=$script_root_dir/plot

mkdir -p $graph_data_dir

t_offset=$(( $t_offset - 1 ))

if [ -f flag_plot_mc ]; then
    echo "folder: $script_plot_dir"

    echo <<H
    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ocn_concat_rg.nc --varname=PDO --colors="$colors" --linestyles="$linestyles" --t-offset=$t_offset --legends=$legends
    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis5_AO.nc --varname=AO --colors="$colors" --linestyles="$linestyles"  --t-offset=$t_offset --legends=$legends
    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ocn_concat_rg.nc --varname=EN34 --colors="$colors" --linestyles="$linestyles"  --t-offset=$t_offset --legends=$legends

H
    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis2.nc --varname=TREFHT_GLB --ylabel="Temperature [K]" --mavg=12 --extra-title="Yearly average" --colors="$colors" --linestyles="$linestyles" --t-offset=$t_offset --legends=$legends
    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis2.nc --varname=TREFHT_LND --ylabel="Temperature [K]" --mavg=12 --extra-title="Yearly average" --colors="$colors" --linestyles="$linestyles" --t-offset=$t_offset --legends=$legends
    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis2.nc --varname=TREFHT_OCN --ylabel="Temperature [K]" --mavg=12 --extra-title="Yearly average" --colors="$colors" --linestyles="$linestyles" --t-offset=$t_offset --legends=$legends
    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ice_trends.nc --varname=total_ice_volume --ylabel="Ice volume [ \$ \\mathrm{km^3} \$ ]" --mavg=12 --yscale="1e9" --colors="$colors" --linestyles="$linestyles" --t-offset=$t_offset --legends=$legends
    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ice_trends.nc --varname=total_ice_area --ylabel="Ice area fraction [ % ]" --mavg=12 --yscale="1e-2" --colors="$colors" --linestyles="$linestyles" --t-offset=$t_offset --legends=$legends

    python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis3_energy.nc --varname=IAET_mean --ylabel="Energy flux [ \$ \\times 10^{15} \\, \\mathrm{W} \$ ]" --yscale="1e15" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles" --legends=$legends
    #python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis4_precip.nc --varname=total_precip_ZM --ylabel="Precipitation [ \$ \\mathrm{mm} \, \\mathrm{yr}^{-1} \$ ]" --yscale="3.171e-11" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends
    #python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis4_precip.nc --varname=total_precip_ZVAR --ylabel="Precipitation [ \$ \\mathrm{mm} \, \\mathrm{yr}^{-1} \$ ]" --yscale="3.171e-11" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends
    
    #python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis4_precip.nc --varname-mean=PREC_TOTAL_ZONAL_MEAN --varname-var=PREC_TOTAL_ZONAL_MAVAR --ylabel="Precipitation [ \$ \\mathrm{mm} \, \\mathrm{yr}^{-1} \$ ]" --yscale="3.171e-11" --domain-file=$atm_domain -colors="$colors" --linestyles="$linestyles"  --legends=$legends

    python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis2.nc --varname=TREFHT_ZM --ylabel="Temperature [ \$ \\mathrm{K} \$ ]" --yscale="1" --y-offset="273.15" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends
    python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis2.nc --varname=TREFHT_ZVAR --ylabel="Temperature [ \$ \\mathrm{K} \$ ]" --yscale="1" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends

    # TREFHT
    python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis2a.nc --varname-mean=TREFHT_ZONAL_MEAN --varname-var=TREFHT_ZONAL_MAVAR --ylabel="Temperature [ \$ \\mathrm{K} \$ ]" --yscale="1" --y-offset="273.15" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="0,:" --extra-title=" Annual mean" --extra-filename="annual_mean" --display-varname="TREFHT" --y-range-mean="-30,30" --y-range-std="0,8"

    # MLD
    python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ocn_analysis3a_rg_MLD.nc --varname-mean=h_ML_ZONAL_MEAN --varname-var=h_ML_ZONAL_MAVAR --ylabel="Mixed-layer Depth [ \$ \\mathrm{m} \$ ]" --yscale="1" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="0,:" --extra-title=" Annual mean" --extra-filename="annual_mean" --display-varname="Mixed-Layer Depth" --y-range-mean="0,400" --y-range-std="0,300"

    # PRECIP    
    python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis4a_precip.nc --varname-mean=PREC_TOTAL_ZONAL_MEAN --varname-var=PREC_TOTAL_ZONAL_MAVAR --ylabel="Precipitation [ \$ \\mathrm{mm} / year \$ ]" --yscale="3.171e-11" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="0,:" --extra-title=" Annual mean" --extra-filename="annual_mean" --display-varname="Precipitation" --y-range-mean="0,4000" --y-range-std="0,1500"

    for m in $( seq 1 12 ); do
        indexing=$( printf "%d,0,:" $(( m - 1 )) )
        python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ocn_analysis3_rg_MLD.nc --varname-mean=h_ML_ZONAL_MEAN --varname-var=h_ML_ZONAL_MAVAR --ylabel="Mixed-layer Depth [ \$ \\mathrm{m} \$ ]" --yscale="1" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="$(( $m - 1 )),0,:" --extra-title=" Month $m" --extra-filename="$(printf '%02d' $m)" --display-varname="Mixed-Layer Depth" --y-range-mean="0,800" --y-range-std="0,300"
        python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis4_precip.nc --varname-mean=PREC_TOTAL_ZONAL_MEAN --varname-var=PREC_TOTAL_ZONAL_MAVAR --ylabel="Precipitation [ \$ \\mathrm{mm} / year \$ ]" --yscale="3.171e-11" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="$(( $m - 1 )),0,:" --extra-title=" Month $m" --extra-filename="$(printf '%02d' $m)" --display-varname="Precipitation" --y-range-mean="0,5000" --y-range-std="0,2000"
    done
fi
