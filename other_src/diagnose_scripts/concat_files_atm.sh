#!/bin/bash

if [ -f $atm_concat ]; then

    echo "$atm_concat already exists. Skip."

else

    echo "Concat atm files of $full_casename"
    cd $atm_hist_dir

    # atm variables
    eval "$(cat <<EOF
    ncrcat -O -v ilev,PSL,V,TREFHT,PRECC,PRECL,FSNT,FSNS,FLNT,FLNS,SHFLX,LHFLX,PSL $full_casename.cam.h0.{$concat_beg_year..$concat_end_year}-{01..12}.nc $atm_concat

EOF
    )"

fi
