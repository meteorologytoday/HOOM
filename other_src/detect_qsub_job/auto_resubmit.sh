#!/bin/bash
#PBS -N piControl_auto_resubmit
#PBS -A UCIR0029
#PBS -l walltime=00:30:00
#PBS -q share
#PBS -j oe
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -m abe


working_dir=/glade/u/home/tienyiao/projects/cesm1_LENS_piControl
status_file=$working_dir/JOB_STATUS

#julia $working_dir/detect_jobs.jl --auto-resubmit --target-year=100 
julia $working_dir/detect_jobs.jl --target-year=74 

job_status=$(tail -n 1 $status_file)
if [ "$job_status" = "UNDONE" ]; then
    echo "Job undone. Sleep for 25 min and resubmit."
    sleep $(( 25 * 60 ))
    qsub $working_dir/auto_resubmit.sh
elif  [ "$job_status" = "DONE" ]; then
    echo "Job all done."
else
    echo "Error: unknown JOB_STATUS: $job_status"
    exit 1
fi
