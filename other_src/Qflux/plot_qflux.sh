#!/bin/bash

export wk_dir=$( dirname $0 )
script_coordtrans_dir=$wk_dir/../CoordTrans



qflux_dir=./Qflux_g16_20191108
transformed_qflux_dir=./Qflux_g16_20191108_transformed
graph_dir=./graph

ocn_domain_file=CESM_domains/domain.ocn.gx1v6.090206.nc
atm_domain_file=CESM_domains/domain.lnd.fv0.9x1.25_gx1v6.090309.nc


# First generate correct transformed  coordinate files
wgt_file=$( basename $ocn_domain_file ".nc" )_$( basename $atm_domain_file ".nc" ).nc

if [ ! -f "$wgt_file" ]; then
    echo "Weight file \"$wgt_file\" does not exist, I am going to generate one..."
    julia -p 4 $script_coordtrans_dir/generate_weight.jl --s-file=$ocn_domain_file --d-file=$atm_domain_file --w-file=$wgt_file --s-mask-value=1.0 --d-mask-value=0.0
fi

for q in $( seq 1 4 ); do

    ocn_qflux_file="$qflux_dir/docn_forcing.Q${q}.energy.g16.nc"
    atm_qflux_file="$transformed_qflux_dir/docn_forcing.Q${q}.energy.f09.nc"
    
    if [ ! -f "$atm_qflux_file" ]; then
        julia $script_coordtrans_dir/transform_data.jl --w-file=$wgt_file --s-file="$ocn_qflux_file" --d-file="$atm_qflux_file" --x-dim=ni --y-dim=nj --t-dim=time --vars=qdp
    fi

    python3 plot_qflux.py --data-file="$atm_qflux_file" --domain-file=$atm_domain_file --output-dir=$graph_dir --output-file-prefix="Q$q" --title-prefix="[Q$q] "

done


