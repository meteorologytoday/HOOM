#!/bin/bash


echo "Transforming data..."

old_files=""
new_files=""

if [ ! -f flag_notrans_ocn ] ; then

    for y in $( seq $diag_beg_year $diag_end_year ); do
        old_file=$full_casename.ocn.h.monthly.$( printf "%04d" $y ).nc
        new_file=$( echo "$old_file" | sed -e 's/\.ocn\./.ocn_rg./' )

        if [ ! -f "$concat_dir/$new_file" ]; then 
            echo "$new_file does not exist, need to transform."
            old_files="${old_file},${old_files}"
            new_files="${new_file},${new_files}"
        fi
    done

    if [ "$old_files" != "" ] ; then

        julia $script_coordtrans_dir/transform_data.jl --s-file=$old_files --d-file=$new_files --w-file=wgt_file.nc --vars=h_ML,T_ML,dTdt_ent,swflx,nswflx,TFLUX_DIV_implied,qflx,qflx2atm,SFLUX_DIV_implied --x-dim=Nx --y-dim=Ny --z-dim=Nz_bone --t-dim=time --s-dir=$ocn_hist_dir --d-dir=$concat_dir
    else
        echo "All files are transformed. Nothing to do here."
    fi
fi


#julia $script_coordtrans_dir/transform_data_ESMF.jl --s-file=$ocn_concat --d-file=$ocn_concat_rg --w-file=wgt_file --vars=h_ML,T_ML,dTdt_ent,swflx,nswflx,TFLUX_DIV_implied,qflx,qflx2atm,SFLUX_DIV_implied --x-dim=Nx --y-dim=Ny --z-dim=Nz_bone --t-dim=time 

