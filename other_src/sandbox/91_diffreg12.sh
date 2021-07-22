#!/bin/bash

./01_run.sh reg1 1 &
./01_run.sh reg5 5 &
wait
ncdiff -O gen_data/record_reg1.nc gen_data/record_reg5.nc regdiff_1_5.nc
