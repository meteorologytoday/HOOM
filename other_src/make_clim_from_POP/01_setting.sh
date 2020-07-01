#!/bin/bash

# skip year one since Jan 01 is not in output
beg_y=2
end_y=21
nyears=$(( end_y - beg_y + 1 ))

wkdir=$(realpath $(dirname $0))

in_dir=$wkdir/ref_run_data
out_dir=$wkdir/data
tmp_dir=$wkdir/tmp
domain_file=$wkdir/CESM_domains/domain.ocn.gx1v6.090206.nc
data_prefix_monthly="paper1_CTL_POP2.pop.h"
data_prefix_daily="paper1_CTL_POP2.pop.h.nday1"
