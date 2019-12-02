#!/bin/bash

echo "Transforming data..."

old_files=""
new_files=""

if [ ! -f flag_notrans_atm ] ; then

    
    for y in $( seq $diag_beg_year $diag_end_year ); do
        for m in $( seq 1 12 ); do
            old_file=$full_casename.cam.h0.$( printf "%04d-%02d" $y $m ).nc
            new_file=$( echo "$old_file" | sed -e 's/\.cam\./.cam_extra./' )

            if [ ! -f "$concat_dir/$new_file" ]; then 
                ((i=i%4)); ((i++==0)) && wait
                echo "$new_file does not exist, need to transform."
                ncap2 -h -O -v -s 'PREC_TOTAL=PRECC+PRECL;' $atm_hist_dir/$old_file $concat_dir/$new_file
            fi
        done
    done

fi

