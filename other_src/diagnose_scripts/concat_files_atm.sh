#!/bin/bash

if [ -f $atm_concat ]; then

    echo "$atm_concat already exists. Skip."

else

    echo "Concat atm files of $full_casename"
    cd $atm_hist_dir

    # atm variables
    eval "$(cat <<EOF
    ncrcat -O -v ilev,PSL,V,TREFHT,PRECC,PRECL,FSNT,FSNS,FLNT,FLNS,SHFLX,LHFLX,PSL $full_casename.cam.h0.{$concat_beg_year..$concat_end_year}-{01..12}.nc $atm_concat
    ncap2 -O -s "PREC_TOTAL=PRECC+PRECL;" $atm_concat $atm_prec
EOF
    )"

fi

if [ -f $atm_prec ]; then

    echo "$atm_prec already exists. Skip."

else

    echo "Creating $atm_prec"
    ncap2 -O -s "PREC_TOTAL=PRECC+PRECL;" $atm_concat $atm_prec

fi
