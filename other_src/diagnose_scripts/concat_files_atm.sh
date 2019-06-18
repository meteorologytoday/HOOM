#!/bin/bash



if [ -f $atm_outputfile ]; then

    echo "No need to concat atm files."
    echo "$atm_outputfile already exists."

else

    printf "Concat atmospheric files... "

    # atm variables
    cd $atm_hist_path 
    eval $(cat <<EOF
    ncrcat -O -v ilev,PSL,V,TREFHT,VQ,VZ,VT,PRECC,PRECL,FSNT,FSNS,FLNT,FLNS,SHFLX,LHFLX,PSL $casename.cam.h0.{$beg_year..$end_year}-{01..12}.nc $atm_outputfile

EOF
    )
    cd $wpath

    printf "done.\n"

fi
