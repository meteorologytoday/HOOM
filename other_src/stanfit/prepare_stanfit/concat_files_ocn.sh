#!/bin/bash

if [ -f $ocn_concat ]; then

    echo "$ocn_concat already exists. Skip."

else

    echo "Concat ocn files of $res_casename"
    cd $ocn_hist_dir

    eval "$(cat <<EOF
    ncrcat -O $res_casename.ocn.h.monthly.{$beg_year..$end_year}.nc $ocn_concat
EOF
    )"

    cd $wdir
    julia $script_coordtrans_dir/transform_data.jl --s-file=$ocn_concat --d-file=$ocn_concat_rg --w-file=$wgt_file --vars=T --x-dim=Nx --y-dim=Ny --z-dim=$Nz --t-dim=time 

fi


