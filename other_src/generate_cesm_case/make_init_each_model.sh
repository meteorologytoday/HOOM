#!/bin/bash

#wk_dir=$( dirname $0 )
script_coordtrans_dir=$wk_dir/../CoordTrans
tmp_dir=tmp

echo "wk_dir: $wk_dir"

lopts=(
    output-file
    label
    data-clim-T-file
    data-clim-S-file
    domain-file
    zdomain-file
    topo-file
    T-unit
    model
    init-config
)

source $wk_dir/getopt_helper.sh

gen_code="$wk_dir/init_code/${model}_${init_config}/make_init.jl"
printf "[%s] => [%s] : [%s]\n" $model $init_config $gen_code

if [ ! -f $output_file ]; then

    julia $wk_dir/init_code/$gen_code               \
        --output-file=$output_file                  \
        --data-clim-T-file=$data_clim_T_file        \
        --data-clim-S-file=$data_clim_S_file        \
        --topo-file=$topo_file                      \
        --domain-file=$domain_file                  \
        --zdomain-file=$zdomain_file                \
        --T-unit=$T_unit

fi
