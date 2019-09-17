#!/bin/bash

lopts=(
    casename
    sim-data-dir
    concat-data-dir
    diag-data-dir
    graph-data-dir
    atm-domain
    ocn-domain
    ice-domain
    concat-beg-year
    concat-end-year
    diag-beg-year
    diag-end-year
    PDO-file
    AO-file
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

export concat_beg_year=$(printf "%04d" $concat_beg_year)
export concat_end_year=$(printf "%04d" $concat_end_year)

export diag_beg_year=$(printf "%04d" $diag_beg_year)
export diag_end_year=$(printf "%04d" $diag_end_year)

export concat_year_stamp=$(printf "%s-%s" $concat_beg_year $concat_end_year)
export diag_year_stamp=$(printf "%s-%s" $diag_beg_year $diag_end_year)

export atm_domain=domain.lnd.fv4x5_gx3v7.091218.nc
export ocn_domain=domain.ocn.gx3v7.120323.nc
export wgt_file=wgt.nc
export PDO_file=PDO_EOFs_fv45.nc
export AO_file=AO_EOFs_fv45.nc


export wdir=`pwd`
export script_root_dir=$(dirname $0)
export script_analysis_dir=$script_root_dir/analysis
export script_coordtrans_dir=$script_root_dir/../CoordTrans
export script_plot_dir=$script_root_dir/plot


echo "Phase 0: Define variables"

export casename=$casename
export full_casename=${casename}

# directories
export concat_dir=$concat_data_dir/$full_casename
export diag_dir=$diag_data_dir/$full_casename
export graph_dir=$graph_data_dir/$full_casename
export sim_dir=$sim_data_dir/$full_casename
export atm_hist_dir=$sim_dir/atm/hist
export ocn_hist_dir=$sim_dir/ocn/hist
export ice_hist_dir=$sim_dir/ice/hist

# filenames
export atm_concat=$concat_dir/atm_concat.nc
export atm_prec=$concat_dir/atm_prec.nc
export ocn_concat=$concat_dir/ocn_concat.nc
export ice_concat=$concat_dir/ice_concat.nc

export ice_concat_rg=$concat_dir/ice_concat_rg.nc
export ocn_concat_rg=$concat_dir/ocn_concat_rg.nc

export atm_analysis1=$diag_dir/atm_analysis1.nc
export atm_analysis2=$diag_dir/atm_analysis2.nc
export atm_analysis3=$diag_dir/atm_analysis3_energy.nc
export atm_analysis4=$diag_dir/atm_analysis4_precip.nc
export atm_analysis5=$diag_dir/atm_analysis5_AO.nc



export ocn_analysis1_rg=$diag_dir/ocn_anomalies1_rg.nc
export ocn_analysis2_rg=$diag_dir/ocn_analysis2_rg_PDO.nc
export ocn_analysis3_rg=$diag_dir/ocn_analysis3_rg_MLD.nc
export ocn_mstat_rg=$diag_dir/ocn_mstat_rg.nc


export ice_analysis1=$diag_dir/ice_trends.nc

for dir_path in $graph_dir $diag_dir $concat_dir ; do

    echo "Checking path: $dir_path"
    mkdir -p $dir_path

done


echo "Phase 1: Generateing weight data, concat and transform data"
begin_t=$SECONDS
if [ ! -f flag_noconcat ]; then

    # Transform ocn grid to atm grid
    if [ ! -f "$wgt_file" ]; then
        echo "Weight file \"$wgt_file\" does not exist, I am going to generate one..."
        julia -p 4  $script_coordtrans_dir/generate_weight.jl --s-file=$ocn_domain --d-file=$atm_domain --w-file=$wgt_file --s-mask-value=1.0 --d-mask-value=0.0
    fi

    echo "Doing case: ${res}_${casename}"
    

    # Doing stuff...
    mkdir -p $diag_dir
    $script_root_dir/concat_files_atm.sh &
    $script_root_dir/concat_files_ocn.sh &
    $script_root_dir/concat_files_ice.sh &

fi
printf "Phase 1 takes %d seconds\n" $(( $SECONDS - $begin_t ))

wait

echo "Phase 2: Doing calculation of diagnose"
begin_t=$SECONDS
echo "Doing case: ${res}_${casename}"

if [ -f flag_diag_all ] || [ -f flag_diag_atm ] ; then
    echo "Diagnose atm..."

    # Need to specify --beg-year --end-year
    julia $script_analysis_dir/atm_anomalies.jl --data-file=$atm_concat --domain-file=$atm_domain --output-file=$atm_analysis1 --beg-year=$diag_beg_year --end-year=$diag_end_year
    julia $script_analysis_dir/atm_temperature.jl --data-file=$atm_concat --domain-file=$atm_domain --output-file=$atm_analysis2 --beg-year=$diag_beg_year --end-year=$diag_end_year
    julia $script_analysis_dir/implied_atm_energy_transport.jl --data-file=$atm_concat --domain-file=$atm_domain --output-file=$atm_analysis3 --beg-year=$diag_beg_year --end-year=$diag_end_year


    julia $script_analysis_dir/mean_anomaly.jl --data-file=$atm_prec --domain-file=$atm_domain --output-file=$atm_analysis4 --beg-year=$diag_beg_year --end-year=$diag_end_year --varname=PREC_TOTAL
   
    echo "Doing averaging over time and lon" 
    #ncwa -O -a month $atm_analysis4 ${atm_analysis4}
    #ncwa -O -a Nx $atm_analysis4 ${atm_analysis4}


    #julia $script_analysis_dir/atm_precip.jl --data-file=$atm_concat --domain-file=$atm_domain --output-file=$atm_analysis4 --beg-year=$diag_beg_year --end-year=$diag_end_year

    # Downstream data. No need to specify --beg-year --end-year
    #julia $script_analysis_dir/AO.jl --data-file=$atm_analysis1 --domain-file=$atm_domain --output-file=$atm_analysis5
fi

if [ -f flag_diag_all ] || [ -f flag_diag_ocn ] ; then

    echo "Diagnose ocn..."
    
    # Need to specify --beg-year --end-year
    #julia $script_analysis_dir/SST_correlation.jl --data-file=$ocn_concat_rg --domain-file=$atm_domain --SST=T_ML --beg-year=$diag_beg_year --end-year=$diag_end_year



    julia $script_analysis_dir/mean_anomaly.jl --data-file=$ocn_concat_rg --domain-file=$atm_domain --output-file=$ocn_analysis3_rg --beg-year=$diag_beg_year --end-year=$diag_end_year --varname=h_ML
    
    # Downstream data. No need to specify --beg-year --end-year
    #julia $script_analysis_dir/PDO.jl --data-file=$ocn_concat_rg --domain-file=$atm_domain --output-file=$ocn_analysis2_rg
    #julia $script_analysis_dir/EN34.jl --data-file-SSTA=$ocn_concat_rg --domain-file=$atm_domain
fi

if [ -f flag_diag_all ] || [ -f flag_diag_ice ] ; then

    echo "Diagnose ice..."
    # Need to specify --beg-year --end-year
    julia $script_analysis_dir/ice.jl --data-file=$ice_concat_rg --domain-file=$atm_domain --output-file=$ice_analysis1 --beg-year=$diag_beg_year --end-year=$diag_end_year

fi

# Plotting
if [ -f flag_plot_sm_all ] || [ -f flag_plot_sm_atm ]; then

    python3 $script_plot_dir/plot_SST.py --data-file=$atm_concat --domain-file=$atm_domain --output-dir=$graph_dir --casename=$casename
    python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casename --data-file=$(basename $atm_analysis3) --varname=IAET_mean --ylabel="Energy flux [ \$ \\times 10^{15} \\, \\mathrm{W} \$ ]" --yscale="1e15" --domain-file=$atm_domain --colors="#000000"
    python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casename --data-file=$(basename $atm_analysis4) --varname=total_precip_ZM --ylabel="Precipitation [ \$ \\mathrm{mm} \, \\mathrm{yr}^{-1} \$ ]" --yscale="3.171e-11" --domain-file=$atm_domain --colors="#000000"
    python3 $script_plot_dir/plot_mc_meridional.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casename --data-file=$(basename $atm_analysis4) --varname=total_precip_ZVAR --ylabel="Precipitation [ \$ \\mathrm{mm} \, \\mathrm{yr}^{-1} \$ ]" --yscale="3.171e-11" --domain-file=$atm_domain --colors="#000000"

    
    python3 $script_plot_dir/plot_contourf.py --data-file=$atm_analysis5 --domain-file=$atm_domain --output-dir=$graph_dir --casename=$casename --varname=PCAs --title="[$casename] PSLA 1st PCA (definition of AO)" --colormap=cmocean_balance --cmin=-1 --cmax=1 --clevs=20 --scale="1" --idx-t=0 --auto-clevs --extra-filename="_AO"
    # Referenced height temperature
    for i in $( seq 1 12 ); do
        echo ""
        python3 $script_plot_dir/plot_contourf.py --data-file=$atm_analysis2 --domain-file=$atm_domain --output-dir=$graph_dir --casename=$casename --varname=TREFHTMM --title="[$casename] Referenced height temperature (month $i)" --colormap=NCV_jet --cmin="-40" --cmax="40" --clevs=20 --idx-t=$(( $i - 1 )) --clabel="[ degC ]" --offset=273.15 --extra-filename="$( printf "_%02d" $i )" --land-transparent
        
    done
fi

if [ -f flag_plot_sm_all ] || [ -f flag_plot_sm_ocn ]; then

    # ocn
    for i in $( seq 1 12 ); do
        echo ""
        python3 $script_plot_dir/plot_ocean_diagnose.py --data-file-SSTAYYC=$ocn_concat_rg --data-file-SSTAVAR=$ocn_concat_rg --domain-file=$atm_domain --output-dir=$graph_dir --casename=$casename --selected-month=$i
    done

    # ocn : PDO
#    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casename --data-file=$(basename $ocn_concat_rg) --varname=PDO --normalize=no --colors="#000000"
    python3 $script_plot_dir/plot_contourf.py --data-file=$ocn_analysis2_rg --domain-file=$atm_domain --output-dir=$graph_dir --casename=$casename --varname=PCAs --title="[$casename] SST 1st PCA (definition of PDO)" --colormap=cmocean_balance --cmin=-1 --cmax=1 --clevs=20 --scale="1" --idx-t=0 --auto-clevs  --extra-filename="_PDO"

fi

if [ -f flag_plot_sm_all ] || [ -f flag_plot_sm_ice ]; then
    # ice
    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casename --data-file=$(basename $ice_analysis1) --varname=total_ice_volume --ylabel="Ice volume [ \$ \\mathrm{km^3} \$ ]" --yscale="1e9" --colors="#000000"
    python3 $script_plot_dir/plot_mc_timeseries.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casename --data-file=$(basename $ice_analysis1) --varname=total_ice_area --ylabel="Ice area fraction [ % ]" --yscale="1e-2" --colors="#000000"

    for i in $( seq 1 12 ); do
        echo ""
        python3 $script_plot_dir/plot_ice.py --data-file=$ice_analysis1 --domain-file=$atm_domain --output-dir=$graph_dir --casename=$casename --selected-month=$i
        python3 $script_plot_dir/plot_contourf.py --data-file=$ice_analysis1 --domain-file=$atm_domain --output-dir=$graph_dir --casename=$casename --varname=aice_MA --title="Sea-ice concentration month $i" --colormap=cmocean_tempo --cmin=10 --cmax=100 --clevs=10 --scale="1e-2" --idx-t=$(( $i - 1 )) --clabel="[ % ]" --extra-filename="$( printf "_%02d" $i )"
    done
fi
printf "Phase 2 takes %d seconds\n" $(( $SECONDS - $begin_t ))


