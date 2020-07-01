#!/bin/bash

beg_y=1
end_y=3
nyears=$(( end_y - beg_y + 1 ))

wkdir=$(realpath $(dirname $0))

in_dir=$wkdir/ref_run_data
out_dir=$wkdir/data
tmp_dir=$wkdir/tmp
domain_file=$wkdir/CESM_domains/domain.ocn.gx1v6.090206.nc
data_prefix_monthly="paper1_CTL_POP2.pop.h"

data_prefix_daily="b.e11.B1850C5CN.f09_g16.005.pop.h.nday1."
data_suffix_daily=".10000101-10991231"

daily_clim_file=$out_dir/daily_clim.nc
