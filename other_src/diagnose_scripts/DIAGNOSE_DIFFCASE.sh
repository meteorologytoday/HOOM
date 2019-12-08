#!/bin/bash
export script_dir=$( dirname "$(realpath $0)" )
export script_plot_dir=$script_dir/plot

# ===== [BEGIN] READ PARAMETERS =====

lopts=(
    case-settings
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


# ===== [END] READ PARAMETERS =====

source "$case_settings"

function join_by { local IFS="$1"; shift; echo "$*"; }

casenames=()
legends=()
colors=()
linestyles=()


ref_casenames=()
ref_legends=()

for i in $(seq 1 $((${#compare_cases[@]}/4))); do
    
    casename=${compare_cases[$((4*(i-1)))]}
    legend=${compare_cases[$((4*(i-1)+1))]}
    color=${compare_cases[$((4*(i-1)+2))]}
    linestyle=${compare_cases[$((4*(i-1)+3))]}

    casenames+=($casename)
    legends+=($legend)
    colors+=($color)
    linestyles+=($linestyle)



    ref_casename=${ref_casenames[$((4*(i-1)))]}
    ref_legend=${ref_casenames[$((4*(i-1)+1))]}
    
    ref_casenames+=($ref_casename)
    ref_legends+=($ref_legend)
    
done

if [ ! -f flag_no_cal_diffcase ] ; then

    for i in $(seq 1 $((${#casenames[@]}))); do
        
        casename=${casenames[$((i-1))]}
        
        printf "Comparing case: [%s] - [%s]\n" $casename $ref_casename
        $script_dir/diff_case.sh --A-dir="${A_dir}/$casename" --B-dir="${B_dir}/${ref_casename}" --AB-dir="${AB_dir}/${casename}_minus_${ref_casename}"
         
    done
fi

if [ -f flag_plot_diffcase ]; then
    mkdir -p $graph_dir    
    for i in $(seq 1 $((${#casenames[@]}))); do
        
        casename=${casenames[$((i-1))]}
        A_casename_dir="${A_dir}/$casename"
        B_casename_dir="${B_dir}/${ref_casename}"
        AB_casename_dir="${AB_dir}/${casename}_minus_${ref_casename}"
        legend=${legends[$((i-1))]}
  
        ref_casename=${ref_casenames[$((i-1))]}
        ref_legend=${ref_legends[$((i-1))]}
 
        # Diff U 
        python3 $script_plot_dir/plot_atm_meridional_mean_std_casediff.py \
            --data-file-1=$A_casename_dir/atm_analysis11_U.nc             \
            --data-file-2=$B_casename_dir/atm_analysis11_U.nc             \
            --domain-file=$atm_domain                                     \
            --output-dir=$graph_dir                                       \
            --varname-mean=U_SM                                           \
            --varname-std=U_SASTD                                         \
            --title="[$legend V.S. $ref_legend] Seasonal analysis of zonal wind" \
            --colormap-mean=bwr                                           \
            --colormap-std=hot_r                                          \
            --colormap-mean-diff=bwr                                      \
            --colormap-std-diff=bwr                                       \
            --cmin-mean=-50                                               \
            --cmax-mean=50                                                \
            --cmax-std=5                                                  \
            --cmax-mean-diff=10                                           \
            --cmax-std-diff=2                                             \
            --clevs-mean=10                                               \
            --clevs-std=5                                                 \
            --clevs-mean-diff=10                                          \
            --clevs-std-diff=4                                            \
            --tick-levs-mean=10                                           \
            --tick-levs-std=5                                             \
            --tick-levs-mean-diff=10                                      \
            --tick-levs-std-diff=4                                        \
            --scale="1"                                                   \
            --extra-filename="U"                                          \
            --idx-x=0                                                     \
            --clabel-mean="Mean [ \$ \\mathrm{m} \, \\mathrm{s}^{-1} \$ ]"\
            --clabel-std="Standard deviation [ \$ \\mathrm{m} \, \\mathrm{s}^{-1} \$ ]"            \
            --clabel-mean-diff="Mean diff [ \$ \\mathrm{m} \, \\mathrm{s}^{-1} \$ ]"               \
            --clabel-std-diff="Standard deviation diff [ \$ \\mathrm{m} \, \\mathrm{s}^{-1} \$ ]"  \
            --lev-file="$A_casename_dir/atm_lev.nc"                                                \
            --subtitles="MAM,JJA,SON,DJF"
    
    
    done
fi


exit;
if [ ! -d "$sim_data_dir" ] ; then
    echo "Error: sim_data_dir='$sim_data_dir' does not exist. "
    exit 1
fi

result_dir=$( printf "%s/result" `pwd`)
concat_data_dir=$( printf "%s/concat" $result_dir )
diag_data_dir=$( printf "%s/%04d-%04d/diag" $result_dir $diag_beg_year $diag_end_year )
graph_data_dir=$( printf "%s/%04d-%04d/graph" $result_dir $diag_beg_year $diag_end_year )

if [ ! -f "$atm_domain" ] ; then
    echo "Error: atm_domain="$atm_domain" does not exists."
    exit 1
fi

if [ ! -f "$ocn_domain" ] ; then
    echo "Error: atm_domain="$ocn_domain" does not exists."
    exit 1
fi


# Parallel loop : https://unix.stackexchange.com/questions/103920/parallelize-a-bash-for-loop

if (( ptasks == 0 )); then
    ptasks=1
fi

echo "### Parallization (ptasks) in batch of $ptasks ###"

# Transform ocn grid to atm grid
if [ ! -f "wgt_file.nc" ]; then
    echo "Weight file \"wgt_file.nc\" does not exist, I am going to generate one..."
    julia -p 4  $script_coordtrans_dir/generate_weight.jl --s-file=$ocn_domain --d-file=$atm_domain --w-file="wgt_file.nc" --s-mask-value=1.0 --d-mask-value=0.0

#    julia $script_coordtrans_dir/generate_SCRIP_format.jl \
#        --input-file=$ocn_domain    \
#        --output-file=SCRIP_${ocn_domain}    \
#        --center-lon=xc     \
#        --center-lat=yc     \
#        --corner-lon=xv     \
#        --corner-lat=yv

#    julia $script_coordtrans_dir/generate_SCRIP_format.jl \
#        --input-file=$atm_domain    \
#        --output-file=SCRIP_${atm_domain}    \
#        --center-lon=xc     \
#        --center-lat=yc     \
#        --corner-lon=xv     \
#        --corner-lat=yv     \
#        --mask-flip

    #ESMF_RegridWeightGen -s SCRIP_${ocn_domain} -d SCRIP_${atm_domain} -m conserve2nd -w $wgt_file --user_areas
#    ESMF_RegridWeightGen -s SCRIP_${ocn_domain} -d SCRIP_${atm_domain} -m neareststod -w $wgt_file --user_areas
fi




for casename in "${casenames[@]}"; do

    ((i=i%ptasks)); ((i++==0)) && wait

    echo "Case: $casename"

    full_casename=${label}_${res}_${casename}
    $script_dir/diagnose_single_model.sh \
        --casename=$casename                \
        --sim-data-dir=$sim_data_dir        \
        --concat-data-dir=$concat_data_dir  \
        --diag-data-dir=$diag_data_dir      \
        --graph-data-dir=$graph_data_dir    \
        --diag-beg-year=$diag_beg_year      \
        --diag-end-year=$diag_end_year      \
        --atm-domain=$atm_domain            \
        --ocn-domain=$ocn_domain            \
        --PCA-sparsity=$PCA_sparsity        & 
done

wait

echo "Start doing model comparison..."

$script_dir/diagnose_mc.sh     \
    --casenames=$( join_by , "${casenames[@]}") \
    --legends=$( join_by , "${legends[@]}") \
    --sim-data-dir=$sim_data_dir                \
    --diag-data-dir=$diag_data_dir              \
    --graph-data-dir=$graph_data_dir            \
    --atm-domain=$atm_domain                    \
    --ocn-domain=$ocn_domain                    \
    --diag-beg-year=$diag_beg_year      \
    --diag-end-year=$diag_end_year      \
    --colors=$( join_by , "${colors[@]}")       \
    --linestyles=$( join_by , "${linestyles[@]}"),     # The comma at the end is necessary. Argparse does not parse "--" as a string however it thinks "--," is a string. 

wait

echo "Done."

