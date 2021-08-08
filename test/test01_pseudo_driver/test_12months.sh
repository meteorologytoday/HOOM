#!/bin/bash

rm -rf Sandbox/archive
mpiexec -n 2 julia --project test22_MPI_HOOM.jl --read-restart=false --stop-n=12 --time-unit=month
ncrcat -O Sandbox/archive/ocn/hist/Sandbox.HOOM.h1.day.0001-{01..12}.nc 1run.nc

rm -rf Sandbox/archive
mpiexec -n 2 julia --project test22_MPI_HOOM.jl --read-restart=false --stop-n=4 --time-unit=month
for i in $( seq 1 2 ); do
    mpiexec -n 2 julia --project test22_MPI_HOOM.jl --read-restart=true --stop-n=4 --time-unit=month
done
ncrcat -O Sandbox/archive/ocn/hist/Sandbox.HOOM.h1.day.0001-{01..12}.nc 4run.nc


ncdiff -O 4run.nc 1run.nc diff.nc
