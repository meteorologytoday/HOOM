#!/bin/bash

if [ -f $ice_concat ] && [ ! -f flag_concat_ice ]; then

    echo "$ice_concat already exists. Skip."

else

    echo "Concat ice files of $full_casename"
    # ice variables
    cd $ice_hist_dir
 
    eval "$(cat <<EOF
    ncrcat -h -O -v aice,hi $full_casename.cice.h.{$concat_beg_year..$concat_end_year}-{01..12}.nc $ice_concat

EOF
    )"

    cd $wdir
    julia $script_coordtrans_dir/transform_data.jl --s-file=$ice_concat --d-file=$ice_concat_rg --w-file=wgt_file.nc --vars=aice,hi --x-dim=ni --y-dim=nj --t-dim=time 

fi
