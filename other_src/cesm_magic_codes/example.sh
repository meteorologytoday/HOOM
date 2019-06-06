#!/bin/bash

. lib_XML.sh

casename=case_example
resolution=f19_g16
mach=cheyenne
compset=E_1850_CN

walltime="12:00:00"
project_code="UCIR0029"

env_run=(
    STOP_OPTION    nyears
    STOP_N         10
    OCN_NCPL       8
    CONTINUE_RUN   FALSE
    RESUBMIT       0
    DOCN_SOM_FILENAME pop_frc.gx3v7.110128.nc
)

env_mach_pes=(
    NTASKS_ATM 70
    ROOTPE_ATM  0
    NTASKS_LND 36
    ROOTPE_LND  0
    NTASKS_ICE 36
    ROOTPE_ICE  0
    NTASKS_OCN  1
    ROOTPE_OCN 71
    NTASKS_CPL 35
    ROOTPE_CPL 36
    NTASKS_GLC  1
    ROOTPE_GLC 69
    NTASKS_ROF 33
    ROOTPE_ROF 36
    NTASKS_WAV  1
    ROOTPE_WAV 70
)

env_vars=(
    caseroot           CASEROOT
    caserun            RUNDIR
    din_loc_root       DIN_LOC_ROOT
    dout_s_root        DOUT_S_ROOT 
    ocn_domain_file    OCN_DOMAIN_FILE
    ocn_domain_path    OCN_DOMAIN_PATH
    totalpes           TOTALPES
    max_tasks_per_node MAX_TASKS_PER_NODE
)


~/ucar_models/cesm1_2_2_1/scripts/create_newcase \
-case      $casename                \
-compset   $compset                     \
-res       $resolution                       \
-mach      $mach


cd $casename



setXML "env_run.xml" "${env_run[@]}"
setXML "env_mach_pes.xml" "${env_mach_pes[@]}"

# Must setup here to get calculated TOTALPES
./cesm_setup

getXML "${env_vars[@]}"

nodes=$(( $totalpes / $max_tasks_per_node ))

# copy user namelist
cp ../../share_nl/user_nl_* .

# Insert config.jl
cat << EOF > config.jl

let

global overwrite_configs = Dict(
    "casename"                   => "$casename",
    "caseroot"                   => "$caseroot",
    "caserun"                    => "$caserun",
    "domain_file"                => "$ocn_domain_path/$ocn_domain_file",
    "short_term_archive_dir"     => "$caserun",
    "long_term_archive_dir"      => "$dout_s_root/ocn",
    "enable_short_term_archive"  => true,
    "enable_long_term_archive"   => true,
    "daily_record"               => true,
    "monthly_record"             => true,
    "yearly_snapshot"            => true,
    "short_term_archive_list"    => "SSM_short_term_archive_list.txt",
    "init_file"                  => "/glade/u/home/tienyiao/projects/cesm1_test/my_inputdata/init_clim_LENS_B1850C5CN_005.nc",
)

end

EOF


cat << EOF > run_ocn.sh

#!/bin/bash
#PBS -N $casename
#PBS -l walltime=$walltime
#PBS -q regular
### Merge output and error files
#PBS -j oe
#PBS -l select=1:ncpus=$max_tasks_per_node:mpiprocs=$max_tasks_per_node
### Send email on abort, begin and end
#PBS -m abe
#PBS -M meteorologytoday@gmail.com

### Run the ocean model ###

LID="\$(date +%y%m%d-%H%M)"
ocn_code="$caseroot/SMARTSLAB-main/src/CESM_driver/run.jl"
config_file="$caseroot/config.jl"

julia -p $max_tasks_per_node \$ocn_code --config="\$config_file" --core=MLMML | tee -a SMARTSLAB.log.\$LID

EOF

cat << EOF > submit_PBS.sh

#!/bin/bash

qsub -A $project_code $caseroot/run_ocn.sh
qsub -A $project_code $caseroot/$casename.run 

EOF

chmod +x run_ocn.sh
chmod +x submit_PBS.sh

# Insert code
git clone https://github.com/meteorologytoday/SMARTSLAB-main.git
cd ./SourceMods/src.docn
ln -s ../../SMARTSLAB-main/src/CESM_driver/cesm1_docn_comp_mod.F90 ./docn_comp_mod.F90
ln -s ../../SMARTSLAB-main/src/CESM_driver/ProgramTunnel .


echo "***** Remember to comment out \`sleep 25\` in $casename.run *****"
echo "***** Remember to change the WALLTIME in $casename.run *****"
