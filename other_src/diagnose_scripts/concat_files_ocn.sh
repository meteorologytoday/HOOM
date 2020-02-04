#!/bin/bash


echo "Transforming data..."



old_files=""
new_files_bilinear=""
new_files_conserve2nd=""

if [ ! -f flag_notrans_ocn ] ; then

    for y in $( seq $diag_beg_year $diag_end_year ); do
        old_file=$full_casename.ocn.h.monthly.$( printf "%04d" $y ).nc
        new_file_bilinear=$( echo "$old_file" | sed -e 's/\.ocn\./.ocn_rg./' )
        new_file_conserve2nd=$( echo "$old_file" | sed -e 's/\.ocn\./.ocn_rg2./' )

        if [ ! -f "$concat_dir/$new_file_bilinear" ] || [ ! -f "$concat_dir/$new_file_conserve2nd" ] ; then 
            echo "$new_file does not exist, need to transform."
            old_files="${old_file},${old_files}"
            new_files_bilinear="${new_file_bilinear},${new_files_bilinear}"
            new_files_conserve2nd="${new_file_conserve2nd},${new_files_conserve2nd}"
        fi
    done

    if [ "$old_files" != "" ] ; then

        #julia $script_coordtrans_dir/transform_data.jl --s-file=$old_files --d-file=$new_files --w-file=wgt_file.nc --vars=h_ML,T_ML,dTdt_ent,swflx,nswflx,TFLUX_DIV_implied,qflx,qflx2atm,SFLUX_DIV_implied,TSAS_clim,SSAS_clim,TFLUX_bot --x-dim=Nx --y-dim=Ny --z-dim=Nz_bone --t-dim=time --s-dir=$ocn_hist_dir --d-dir=$concat_dir
        julia $script_coordtrans_dir/transform_data.jl --s-file=$old_files --d-file=$new_files_bilinear --w-file=${wgt_dir}/wgt.bilinear.nc --algo=ESMF --vars=h_ML,T_ML --x-dim=Nx --y-dim=Ny --z-dim=Nz_bone --t-dim=time --s-dir=$ocn_hist_dir --d-dir=$concat_dir
        julia $script_coordtrans_dir/transform_data.jl --s-file=$old_files --d-file=$new_files_conserve2nd --w-file=${wgt_dir}/wgt.conserve2nd.nc --algo=ESMF --vars=TFLUX_DIV_implied --x-dim=Nx --y-dim=Ny --t-dim=time --s-dir=$ocn_hist_dir --d-dir=$concat_dir
    else
        echo "All files are transformed. Nothing to do here."
    fi
fi


#julia $script_coordtrans_dir/transform_data_ESMF.jl --s-file=$ocn_concat --d-file=$ocn_concat_rg --w-file=wgt_file --vars=h_ML,T_ML,dTdt_ent,swflx,nswflx,TFLUX_DIV_implied,qflx,qflx2atm,SFLUX_DIV_implied --x-dim=Nx --y-dim=Ny --z-dim=Nz_bone --t-dim=time 

