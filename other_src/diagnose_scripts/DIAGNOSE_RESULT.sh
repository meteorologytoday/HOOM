#!/bin/bash


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


# Output kill process shell.
echo "$(cat <<EOF
#!/bin/bash
kill -9 -$$
EOF
)" > kill_process.sh
chmod +x kill_process.sh

function join_by { local IFS="$1"; shift; echo "$*"; }



casenames=()
legends=()
colors=()
linestyles=()

for i in $(seq 1 $((${#case_settings[@]}/4))); do
    casename=${case_settings[$((4*(i-1)))]}
    legend=${case_settings[$((4*(i-1)+1))]}
    color=${case_settings[$((4*(i-1)+2))]}
    linestyle=${case_settings[$((4*(i-1)+3))]}
    printf "[%s] => [%s, %s]\n" $casename $color $linestyle

    casenames+=($casename)
    legends+=($legend)
    colors+=($color)
    linestyles+=($linestyle)
done


if [ ! -d "$sim_data_dir" ] ; then
    echo "Error: sim_data_dir='$sim_data_dir' does not exist. "
    exit 1
fi

result_dir=$( printf "%s/result_%04d-%04d" `pwd` $concat_beg_year $concat_end_year )
concat_data_dir=$( printf "%s/concat" $result_dir )
diag_data_dir=$( printf "%s/%04d-%04d/diag" $result_dir $diag_beg_year $diag_end_year )
graph_data_dir=$( printf "%s/%04d-%04d/graph" $result_dir $diag_beg_year $diag_end_year )

atm_domain=domain.lnd.fv4x5_gx3v7.091218.nc
ocn_domain=domain.ocn.gx3v7.120323.nc



# Parallel loop : https://unix.stackexchange.com/questions/103920/parallelize-a-bash-for-loop

if (( ptasks == 0 )); then
    ptasks=1
fi

echo "### Parallization (ptasks) in batch of $ptasks ###"

for casename in "${casenames[@]}"; do

    ((i=i%N)); ((i++==0)) && wait

    echo "Case: $casename"

    full_casename=${label}_${res}_${casename}
     ./other_src/diagnose_scripts/diagnose_single_model.sh \
        --casename=$casename                \
        --sim-data-dir=$sim_data_dir        \
        --concat-data-dir=$concat_data_dir  \
        --diag-data-dir=$diag_data_dir      \
        --graph-data-dir=$graph_data_dir    \
        --concat-beg-year=$concat_beg_year  \
        --concat-end-year=$concat_end_year  \
        --diag-beg-year=$diag_beg_year      \
        --diag-end-year=$diag_end_year      \
        --atm-domain=$atm_domain            \
        --ocn-domain=$ocn_domain            & 
done

wait

echo "Start doing model comparison..."

./other_src/diagnose_scripts/diagnose_mc.sh     \
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

