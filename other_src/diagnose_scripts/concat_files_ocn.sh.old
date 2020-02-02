#!/bin/bash

if [ -f $ocn_concat_rg ] && [ ! -f flag_concat_ocn ]; then

    echo "$ocn_concat_rg already exists. Skip."

else

    if [ -f flag_concat_ocn ] || [ ! -f $ocn_concat ]; then

        echo "Concat ocn files of $full_casename"
        cd $ocn_hist_dir

        eval "$(cat <<EOF
        ncrcat -h -O -v h_ML,T_ML,dTdt_ent,swflx,nswflx,TFLUX_DIV_implied,qflx,qflx2atm,SFLUX_DIV_implied $full_casename.ocn.h.monthly.{$concat_beg_year..$concat_end_year}.nc $ocn_concat
EOF
        )"
        #ncrcat -O -v h_ML,T_ML,dTdt_ent,T_hadvs,T_vadvs,swflx,nswflx $full_casename.ocn.h.monthly.{$concat_beg_year..$concat_end_year}.nc $ocn_concat
        #rm tmp.ocn.*.nc
        #for y in $(seq $concat_beg_year $concat_end_year); do
        #    echo "YEAR $y"
        #    ncks -h -O -d Nz_bone,0,0 -v T_hadvs,T_vadvs "$full_casename.ocn.h.monthly.$(printf '%04d' $y).nc" tmp.ocn.1.nc
        #    if [ -f tmp.ocn.2.nc ]; then
        #        ncrcat -h --rec_apn tmp.ocn.1.nc tmp.ocn.2.nc
        #    else
        #        mv tmp.ocn.1.nc tmp.ocn.2.nc
        #    fi
        #done
        #ncks -h -A tmp.ocn.2.nc $ocn_concat

        #rm tmp.ocn.*.nc        
    fi

    echo "Transforming data..."

    cd $wdir
    #julia $script_coordtrans_dir/transform_data_ESMF.jl --s-file=$ocn_concat --d-file=$ocn_concat_rg --w-file=wgt_file --vars=h_ML,T_ML,dTdt_ent,swflx,nswflx,TFLUX_DIV_implied,qflx,qflx2atm,SFLUX_DIV_implied --x-dim=Nx --y-dim=Ny --z-dim=Nz_bone --t-dim=time 
    julia $script_coordtrans_dir/transform_data.jl --s-file=$ocn_concat --d-file=$ocn_concat_rg --w-file=wgt_file.nc --vars=h_ML,T_ML,dTdt_ent,swflx,nswflx,TFLUX_DIV_implied,qflx,qflx2atm,SFLUX_DIV_implied --x-dim=Nx --y-dim=Ny --z-dim=Nz_bone --t-dim=time 

fi


