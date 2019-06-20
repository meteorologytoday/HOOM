#!/bin/bash

lopts=(res casename sim-data-dir diag-data-dir graph-dir atm-domain ocn-domain beg-year end-year  PDO-file AO-file)

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

export beg_year=$(printf "%04d" $beg_year)
export end_year=$(printf "%04d" $end_year)
export year_stamp=$(printf "%s-%s" $beg_year $end_year)

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

for dir_path in $graph_dir $diag_data_dir ; do

    echo "Checking path: $dir_path"
    mkdir -p $dir_path

done

echo "Phase 0: Define variables"

if [[ $casename =~ SSM_NK ]]; then
    export Nz=Nz_bone
else
    export Nz=Nz
fi

export casename=$casename
export res_casename=${res}_${casename}

# directories
export diag_dir=$diag_data_dir/$res_casename
export sim_dir=$sim_data_dir/$res_casename
export atm_hist_dir=$sim_dir/atm/hist
export ocn_hist_dir=$sim_dir/ocn

# filenames
export atm_concat=$diag_dir/atm_concat.nc
export ocn_concat=$diag_dir/ocn_concat.nc

export atm_analysis1=$diag_dir/atm_analysis1.nc
export atm_analysis2=$diag_dir/atm_analysis2.nc

export ocn_concat_rg=$diag_dir/ocn_concat_rg.nc

export ocn_analysis1_rg=$diag_dir/ocn_anomalies1_rg.nc
export ocn_mstat_rg=$diag_dir/ocn_mstat_rg.nc




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

fi
printf "Phase 1 takes %d seconds\n" $(( $SECONDS - $begin_t ))

wait

echo "Phase 2: Doing calculation of diagnose"
begin_t=$SECONDS
echo "Doing case: ${res}_${casename}"

if [ ! -f flag_nodiag ]; then
    $script_root_dir/diagnose_calculation.sh
fi

if [ ! -f flag_noplot ]; then

    # atm
    python3 $script_plot_dir/plot_SST.py --data-file=$atm_concat --domain-file=$atm_domain --output-dir=$graph_dir --casename=$casename

    # ocn
    for i in $( seq 1 12 ); do
        python3 $script_plot_dir/plot_ocean_diagnose.py --data-file-SSTAYYC=$ocn_concat_rg --data-file-SSTAVAR=$ocn_concat_rg --domain-file=$atm_domain --output-dir=$graph_dir --casename=$casename --selected-month=$i
    done

    # PDO
    python3 $script_plot_dir/plot_mc_climate_indices.py --input-dir=$diag_data_dir --output-dir=$graph_dir --res=$res --casenames=$casename --data-file=$(basename $ocn_concat_rg) --varname=PDO --normalize=no

fi
printf "Phase 2 takes %d seconds\n" $(( $SECONDS - $begin_t ))

