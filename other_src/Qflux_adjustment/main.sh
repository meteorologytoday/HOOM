#!/bin/bash
lopts=(
    model
    fix
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



casename="LENS_Qflux_finding_with_IFRAC_f09_g16_${model}_20"
zero_qflx="../raw_data/docn_forcing.20200601.g16.daily.init.ConH.nc"
qflx_fix_fileformat="../raw_data/docn_forcing.g16.daily.%s.fixed%02d.nc"

year_unit=10
beg_year=$(( ( fix - 1 ) * year_unit + 1 ))
end_year=$(( beg_year + year_unit - 1 ))

data_file_prefix="/glade/u/home/tienyiao/scratch-tienyiao/archive/${casename}/ocn/hist/${casename}.ocn.h.daily."
domain_file="../raw_data/domain.ocn.gx1v6.090206.nc"
correction_file=$( printf "../raw_data/qflx_correction_%s_%02d-%02d.nc" $model $beg_year $end_year )
new_qflx_file=$( printf "$qflx_fix_fileformat" $model $fix )

if [ "$fix" -eq "1" ]; then
    old_qflx_file="$zero_qflx"
else
    old_qflx_file=$( printf "$qflx_fix_fileformat" $model $(( fix - 1 )) )
fi

./do_correction_daily.sh  \
    --beg-year=$beg_year                   \
    --end-year=$end_year                   \
    --domain-file=$domain_file             \
    --data-file-prefix="$data_file_prefix" \
    --old-qflx-file="$old_qflx_file"       \
    --new-qflx-file="$new_qflx_file"       \
    --correction-file="$correction_file" 


