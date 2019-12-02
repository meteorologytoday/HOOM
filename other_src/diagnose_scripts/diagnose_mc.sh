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

    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ocn_analysis2_rg_PDO.nc --varname=PDO --colors="$colors" --linestyles="$linestyles" --t-offset=$t_offset --legends=$legends --mavg=6
    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis5_AO.nc --varname=AO --colors="$colors" --linestyles="$linestyles"  --t-offset=$t_offset --legends=$legends --mavg=6
    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ocn_analysis11_rg_EN34.nc --varname=EN34 --colors="$colors" --linestyles="$linestyles"  --t-offset=$t_offset --legends=$legends --mavg=6
    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ocn_analysis11_rg_EN34.nc --varname=ENSO --colors="$colors" --linestyles="$linestyles"  --t-offset=$t_offset --legends=$legends --mavg=6

    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis7.nc --varname=TREFHT_GLB --ylabel="Temperature [deg C]" --mavg=12 --extra-title="Yearly average" --colors="$colors" --linestyles="$linestyles" --t-offset=$t_offset --y-offset="273.15" --legends=$legends
    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis7.nc --varname=TREFHT_LND --ylabel="Temperature [deg C]" --mavg=12 --extra-title="Yearly average" --colors="$colors" --linestyles="$linestyles" --t-offset=$t_offset --y-offset="273.15" --legends=$legends
    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis7.nc --varname=TREFHT_OCN --ylabel="Temperature [deg C]" --mavg=12 --extra-title="Yearly average" --colors="$colors" --linestyles="$linestyles" --t-offset=$t_offset --y-offset="273.15" --legends=$legends


    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ice_trends.nc --varname=total_ice_volume --ylabel="Ice volume [ \$ \\mathrm{km^3} \$ ]" --mavg=12 --yscale="1e9" --colors="$colors" --linestyles="$linestyles" --t-offset=$t_offset --legends=$legends
    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ice_trends.nc --varname=total_ice_area --ylabel="Ice area fraction [ % ]" --mavg=12 --yscale="1e-2" --colors="$colors" --linestyles="$linestyles" --t-offset=$t_offset --legends=$legends

    #python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis3_energy.nc --varname=IAET_mean --ylabel="Energy flux [ \$ \\times 10^{15} \\, \\mathrm{W} \$ ]" --yscale="1e15" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles" --legends=$legends
    
    python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis3_energy.nc --varname=IAET_mean --ylabel="Energy flux [ \$ \\times 10^{15} \\, \\mathrm{W} \$ ]" --yscale="1e15" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles" --legends=$legends
#    python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ocn_analysis9b_rg_oiet_avg.nc --varname=IET_OCN_ZONAL_MEAN --ylabel="Energy flux [ \$ \\times 10^{15} \\, \\mathrm{W} \$ ]" --yscale="1e15" --domain-file=$atm_domain --colors="$colors" --indexing="0,:" --linestyles="$linestyles" --legends=$legends

#    python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis3_energy.nc --varname=L_IET --ylabel="Energy flux [ \$ \\times 10^{15} \\, \\mathrm{W} \$ ]" --yscale="1e15" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles" --legends=$legends

    #python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis4_precip.nc --varname=total_precip_ZM --ylabel="Precipitation [ \$ \\mathrm{mm} \, \\mathrm{yr}^{-1} \$ ]" --yscale="3.171e-11" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends
    #python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis4_precip.nc --varname=total_precip_ZVAR --ylabel="Precipitation [ \$ \\mathrm{mm} \, \\mathrm{yr}^{-1} \$ ]" --yscale="3.171e-11" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends
    
    #python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis4_precip.nc --varname-mean=PREC_TOTAL_ZONAL_MEAN --varname-var=PREC_TOTAL_ZONAL_MAVAR --ylabel="Precipitation [ \$ \\mathrm{mm} \, \\mathrm{yr}^{-1} \$ ]" --yscale="3.171e-11" --domain-file=$atm_domain -colors="$colors" --linestyles="$linestyles"  --legends=$legends

    python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis2.nc --varname=TREFHT_ZM --ylabel="Temperature [ \$ \\mathrm{K} \$ ]" --yscale="1" --y-offset="273.15" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends
    python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis2.nc --varname=TREFHT_ZVAR --ylabel="Temperature [ \$ \\mathrm{K} \$ ]" --yscale="1" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends

    # TREFHT
    python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis2a.nc --varname-mean=TREFHT_ZONAL_MEAN --varname-var=TREFHT_ZONAL_MAVAR --ylabel="Temperature [ \$ \\mathrm{K} \$ ]" --yscale="1" --y-offset="273.15" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="0,:" --extra-title=" Annual mean" --extra-filename="annual_mean" --display-varname="TREFHT" --y-range-mean="-30,30" --y-range-std="0,8"

    # MLD
   python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ocn_analysis3a_rg_MLD.nc --varname-mean=h_ML_ZONAL_MEAN --varname-var=h_ML_ZONAL_MAVAR --ylabel="Mixed-layer Depth [ \$ \\mathrm{m} \$ ]" --yscale="1" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="0,:" --extra-title=" Annual mean" --extra-filename="annual_mean" --display-varname="Mixed-Layer Depth" --y-range-mean="0,400" --y-range-std="0,300" --invert-y-axis

    # PRECIP    
    python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis4a_precip.nc --varname-mean=PREC_TOTAL_ZONAL_MEAN --varname-var=PREC_TOTAL_ZONAL_MAVAR --ylabel="Precipitation [ \$ \\mathrm{mm} / year \$ ]" --yscale="3.171e-11" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="0,:" --extra-title=" Annual mean" --extra-filename="annual_mean" --display-varname="Precipitation" --y-range-mean="0,4000" --y-range-std="0,1500"

    # ICEFRAC
    python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis6a_icefrac.nc --varname-mean=ICEFRAC_ZONAL_MEAN --varname-var=ICEFRAC_ZONAL_MAVAR --ylabel="Concentration [ \$ \\% \$ ]" --yscale="1e-2" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="0,:" --extra-title=" Annual mean" --extra-filename="annual_mean" --display-varname="Sea-ice Fraction" --y-range-mean="0,100" --y-range-std="0,50"

    # Solar surface flux
    python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ocn_analysis7a_rg_swflx.nc --varname-mean=swflx_ZONAL_MEAN --varname-var=swflx_ZONAL_MAVAR --ylabel="Energy density flux [ \$ \\mathrm{W} \\, \\mathrm{m}^{-2} \$ ]" --yscale="-1.0" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="0,:" --extra-title=" Annual mean" --extra-filename="annual_mean" --display-varname="Solar surface flux" --y-range-mean="0,300" --y-range-std="0,50"

    # Non-solar surface flux
    python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ocn_analysis8a_rg_nswflx.nc --varname-mean=nswflx_ZONAL_MEAN --varname-var=nswflx_ZONAL_MAVAR --ylabel="Energy density flux [ \$ \\mathrm{W} \\, \\mathrm{m}^{-2} \$ ]" --yscale="-1.0" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="0,:" --extra-title=" Annual mean" --extra-filename="annual_mean" --display-varname="Non-solar surface flux" --y-range-mean="-200, 100" --y-range-std="0,50"

    # Ekman dTdt
    python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ocn_analysis5a_rg_T_hadvs.nc --varname-mean=T_hadvs_ZONAL_MEAN --varname-var=T_hadvs_ZONAL_MAVAR --ylabel="Temperature tendency [ \$ \\mathrm{K} \\, \\mathrm{mon}^{-1} \$ ]" --yscale="3.858e-7" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="0,:" --extra-title=" Annual mean" --extra-filename="annual_mean" --display-varname="Temperature change caused by horizontal advection" --y-range-mean="-1, 1" --y-range-std="0,1"



    for m in $( seq 1 12 ); do
        indexing=$( printf "%d,0,:" $(( m - 1 )) )

    

        # TREFHT
        python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis2.nc --varname-mean=TREFHT_ZONAL_MEAN --varname-var=TREFHT_ZONAL_MAVAR --ylabel="Temperature [ \$ \\mathrm{deg C} \$ ]" --yscale="1" --y-offset="273.15" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="$(( $m - 1 )),0,:" --extra-title=" Month $m" --extra-filename="$(printf '%02d' $m)" --display-varname="Referenced Height Temperature" --y-range-mean="-30,30" --y-range-std="0,8"


        # h_ML
        python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ocn_analysis3_rg_MLD.nc --varname-mean=h_ML_ZONAL_MEAN --varname-var=h_ML_ZONAL_MAVAR --ylabel="Mixed-layer Depth [ \$ \\mathrm{m} \$ ]" --yscale="1" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="$(( $m - 1 )),0,:" --extra-title=" Month $m" --extra-filename="$(printf '%02d' $m)" --display-varname="Mixed-Layer Depth" --y-range-mean="0,800" --y-range-std="0,300" --invert-y-axis

        # PRECIP
        python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis4_precip.nc --varname-mean=PREC_TOTAL_ZONAL_MEAN --varname-var=PREC_TOTAL_ZONAL_MAVAR --ylabel="Precipitation [ \$ \\mathrm{mm} / year \$ ]" --yscale="3.171e-11" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="$(( $m - 1 )),0,:" --extra-title=" Month $m" --extra-filename="$(printf '%02d' $m)" --display-varname="Precipitation" --y-range-mean="0,5000" --y-range-std="0,2000"
        
        # ICEFRAC
        python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=atm_analysis6_icefrac.nc --varname-mean=ICEFRAC_ZONAL_MEAN --varname-var=ICEFRAC_ZONAL_MAVAR --ylabel="Concentration [ \$ \\% \$ ]" --yscale="1e-2" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="$(( $m - 1)),0,:" --extra-title=" Month $m" --extra-filename="$(printf '%02d' $m)" --display-varname="Sea-ice Fraction" --y-range-mean="0,100" --y-range-std="0,50"

        # dTdt horizontal advection
        python3 $script_plot_dir/plot_mc_meridional_mean_std.py --input-dir=$diag_data_dir --output-dir=$graph_data_dir --casenames=$casenames --data-file=ocn_analysis5_rg_T_hadvs.nc --varname-mean=T_hadvs_ZONAL_MEAN --varname-var=T_hadvs_ZONAL_MAVAR --ylabel="Temperature Tendency [ \$ \\mathrm{K} \\, \\mathrm{mon}^{-1} \$ ]" --yscale="3.858e-7" --domain-file=$atm_domain --colors="$colors" --linestyles="$linestyles"  --legends=$legends --indexing="$(( $m - 1)),0,:" --extra-title=" Month $m" --extra-filename="$(printf '%02d' $m)" --display-varname="Temperature Tendency caused by horizontal advection" --y-range-mean="-1, 1" --y-range-std="0,1"
 

    done
fi
