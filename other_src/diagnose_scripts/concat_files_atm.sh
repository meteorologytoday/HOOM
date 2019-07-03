#!/bin/bash

if [ -f $atm_concat ]; then

    echo "$atm_concat already exists. Skip."

else

    echo "Concat atm files of $res_casename"

    # atm variables
    cd $atm_hist_dir 
    eval "$(cat <<EOF
    ncrcat -O -v ilev,PSL,V,TREFHT,VQ,VZ,VT,PRECC,PRECL,FSNT,FSNS,FLNT,FLNS,SHFLX,LHFLX,PSL $res_casename.cam.h0.{$beg_year..$end_year}-{01..12}.nc $atm_concat

EOF
    )"

fi
