#!/bin/bash

export wk_dir=$( dirname $0 )
script_coordtrans_dir=$wk_dir/../CoordTrans



graph_dir=./graph

ocn_domain_file=CESM_domains/domain.ocn.gx1v6.090206.nc
atm_domain_file=CESM_domains/domain.lnd.fv0.9x1.25_gx1v6.090309.nc


# First generate correct transformed  coordinate files
wgt_file=$( basename $ocn_domain_file ".nc" )_$( basename $atm_domain_file ".nc" ).nc

if [ ! -f "$wgt_file" ]; then
    echo "Weight file \"$wgt_file\" does not exist, I am going to generate one..."
    julia -p 4 $script_coordtrans_dir/generate_weight.jl --s-file=$ocn_domain_file --d-file=$atm_domain_file --w-file=$wgt_file --s-mask-value=1.0 --d-mask-value=0.0
fi


mkdir -p $graph_dir


ocn_seaice_file="dice_forcing.gx1v6.piControl.nc"
atm_seaice_file="dice_forcing.f09.piControl.nc"
 
#ocn_seaice_file="dice_forcing.gx1v6.icefree.nc"
#atm_seaice_file="dice_forcing.f09.icefree.nc"
    
if [ ! -f "$atm_seaice_file" ]; then
    julia $script_coordtrans_dir/transform_data.jl --w-file=$wgt_file --s-file="$ocn_seaice_file" --d-file="$atm_seaice_file" --x-dim=lon --y-dim=lat --t-dim=time --vars=ice_cov
fi

python3 plot_seaice.py --data-file="$atm_seaice_file" --domain-file=$atm_domain_file --output-dir=$graph_dir --output-file-prefix="piControl" --title-prefix="[CTL] "
#python3 plot_seaice.py --data-file="$atm_seaice_file" --domain-file=$atm_domain_file --output-dir=$graph_dir --output-file-prefix="icemelt" --title-prefix="[ICEMELT] "

