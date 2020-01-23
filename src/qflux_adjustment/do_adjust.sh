#!/bin/bash

julia qflux_adjustment.jl \
    --input-file=data/docn_forcing.Q1.energy.g16.nc   \
    --output-file=data/updated_qflux.nc               \
    --target-SST-file=data/SST_monthly_clim.nc             \
    --current-SST-file=data/SST_SOM.nc   
