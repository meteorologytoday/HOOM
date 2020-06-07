#!/bin/bash

lopts=(
    beg-year
    end-year
    domain-file
    data-file-prefix
    old-qflx-file
    new-qflx-file
    correction-file
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

export script_dir=$( dirname "$(realpath $0)" )

julia $script_dir/mk_correction_file_daily.jl        \
    --beg-year=$beg_year                   \
    --end-year=$end_year                   \
    --domain-file=$domain_file             \
    --data-file-prefix="$data_file_prefix" \
    --input-file="$old_qflx_file"          \
    --output-file="$correction_file"       

cp $old_qflx_file tmp.nc
ncks -A -v Qflx_T_final,Qflx_S_final $correction_file tmp.nc
ncap2 -O -s 'Qflx_T=Qflx_T_final;Qflx_S=Qflx_S_final;' tmp.nc tmp.nc
ncks -O -x -v Qflx_T_final,Qflx_S_final tmp.nc $new_qflx_file

rm tmp.nc

