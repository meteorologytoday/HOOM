# SMARTSLAB

This project aims to develope mixed-layer based ocean component. It also extends the data ocean model (docn) of Community Earth System Model (CESM). It is an experimental framework for fortran based program to work with Julia program without deep embedding.


# Usage

## CESM 1.2.2.1

### Configuration
After creating a CESM case, this is what typically needs to be done.

1. Copy or create soft link of `docn_comp_mod.F90` and `fortran_lib` to `$CASEROOT/SourceMods/src.docn`. The content of `docn_comp_mod.F90` should be the file `cesm1_comp_mod.F90`.
2. Edit `$CASEROOT/env_run.xml` entries:
    - STOP_N : xxx (Days of simulation you want)
    - OCN_NCPL : 1 (Ocean update frequency: 1 time per day )
    - DOCN_SOM_FILENAME : pop_frc.gx3v7.110128.nc (actually any valid q-flux file but remember the data will not be used. This is only to make docn run).
3. Edit `$CASEROOT/env_mach_pes.xml` entries (Please study [CESM1 PES LAYOUT GUIDE](http://www.cesm.ucar.edu/models/cesm1.0/cesm/cesm_doc_1_0_6/x2565.html)) :
    - MAX_TASKS_PER_NODE : 4
    - NTASKS_OCN : 1 (Must be 1)
    - NTHRDS_OCN : 1
4. Create your own `config.jl`(Here after `$SMARTSLAB_CONFIG`) file for ocean model. You need to specify `caseroot`, `wdir` ($CASERUN), `domain_file` (ex: `$CESMINPUTDATA/share/domains/domain.ocn.gx3v7.120323.nc`)
5. Run `$CASEROOT/cesm_setup`, then `$CASEROOT/$CASE.build`.
6. Edit `$CASEROOT/$CASE.run`:
    - Remove the `sleep 25` command (ill!)
    - Uncomment the line of either `mpirun` or `mpiexec` depending on which command is available to you.
7. Create named-pipe files: `$CASERUN/cesm2mymodel.pipe` and `$CASERUN/mymodel2cesm.pipe`.

### Run
1. Run `$CASEROOT/$CASE.run`
2. Run `julia $SMARTSLABDEV_ROOT/src/CESM_driver/CESM_MLMML/run_CESM_MLMML.jl --config=`$SMARTSLAB_CONFIG`

