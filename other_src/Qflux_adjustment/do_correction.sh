#!/bin/bash

lopts=(
    old-qflux-file
    correction-file
    new-qflux-file
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

cp $old_qflux_file tmp.nc
ncks -A qflx_correction_mean $correction_file tmp.nc
ncap2 -O -s 'qdp = qdp + qflx_correction_mean' tmp.nc tmp.nc
ncks -x -v qflx_correction_mean $new_qflux_file

