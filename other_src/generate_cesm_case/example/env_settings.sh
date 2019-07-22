#!/bin/bash

env_run=(
    STOP_OPTION       nmonths
    STOP_N            1
    OCN_NCPL          1 
    CONTINUE_RUN      FALSE
    RESUBMIT          1
    DOCN_FILENAME     "pop_frc.gx3v7.110128.nc"
)

env_mach_pes=(
    NTASKS_ATM 35
    ROOTPE_ATM  0
    NTASKS_LND 35
    ROOTPE_LND  0
    NTASKS_ICE 35
    ROOTPE_ICE  0
    NTASKS_OCN  1
    ROOTPE_OCN 35
    NTASKS_CPL 35
    ROOTPE_CPL  0
    NTASKS_GLC 35
    ROOTPE_GLC  0
    NTASKS_ROF 35
    ROOTPE_ROF  0
    NTASKS_WAV 35
    ROOTPE_WAV  0
)

env_mach_pes=(
    NTASKS_ATM  3
    ROOTPE_ATM  0
    NTASKS_LND  3
    ROOTPE_LND  0
    NTASKS_ICE  3
    ROOTPE_ICE  0
    NTASKS_OCN  1
    ROOTPE_OCN  3
    NTASKS_CPL  3
    ROOTPE_CPL  0
    NTASKS_GLC  3
    ROOTPE_GLC  0
    NTASKS_ROF  3
    ROOTPE_ROF  0
    NTASKS_WAV  3
    ROOTPE_WAV  0
)


