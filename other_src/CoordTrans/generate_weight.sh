#!/bin/bash

script_root_dir=$(dirname $0)

lopts=(
    s-file
    d-file
    w-file
    s-mask-value
    d-mask-value
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

s_SCRIP=SCRIP_$( basename ${s_file} )
d_SCRIP=SCRIP_$( basename ${d_file} )

julia $script_root_dir/generate_SCRIP_format.jl \
    --input-file=${s_file}    \
    --output-file=${s_SCRIP}    \
    --center-lon=xc     \
    --center-lat=yc     \
    --corner-lon=xv     \
    --corner-lat=yv     \
    --mask-value=${s_mask_value}

julia $script_root_dir/generate_SCRIP_format.jl \
    --input-file=${d_file}    \
    --output-file=${d_SCRIP}    \
    --center-lon=xc     \
    --center-lat=yc     \
    --corner-lon=xv     \
    --corner-lat=yv     \
    --mask-value=${d_mask_value}

ESMF_RegridWeightGen -s ${s_SCRIP} -d ${d_SCRIP} -m conserve2nd -w $w_file --user_areas --check




