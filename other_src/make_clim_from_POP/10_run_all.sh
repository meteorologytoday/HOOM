#!/bin/bash



./02_make_clim_mean_and_init.sh
./03_cal_daily_from_monthly_SSS_HMXL.sh
./04_cal_daily_mean_SST.sh
./05_generate_daily_SOM_qflx.sh
./06_generate_final_file.sh

