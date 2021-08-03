#!/bin/bash

mpiexec -n 2 julia --project test22_MPI_HOOM.jl
cp Sandbox/caserun/Sandbox.HOOM.h1.day.0001-01.nc result_core02.nc


mpiexec -n 5 julia --project test22_MPI_HOOM.jl
cp Sandbox/caserun/Sandbox.HOOM.h1.day.0001-01.nc result_core05.nc

ncdiff -O result_core05.nc  result_core02.nc diff_0502.nc
