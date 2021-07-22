#!/bin/bash

source 00_path.sh

set -x
julia -p ${procs} $wkdir/run_model.jl \
    --init-file=$init_file            \
    --output-file=$record_file        \
    --domain-file=$domain_file        \
    --run-days=50

