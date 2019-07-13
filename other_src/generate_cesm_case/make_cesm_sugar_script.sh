#!/bin/bash

wk_dir=$( dirname $0 )
echo $( basename $0 )
echo "wk_dir: $wk_dir"

lopts=(
    code-output-dir
    init-file
    label
    resolution
    walltime
    project-code
    compset
    machine
    cesm-create-newcase
    user-namelist-dir
    model
)

source $wk_dir/getopt_helper.sh


export casename="${resolution}_${label}"

script_file=$code_output_dir/make_case_$casename.sh


cat $wk_dir/lib_XML.sh >> $script_file

cat << EOF >> $script_file

casename=$casename
resolution=$resolution
machine=$machine
compset=$compset
user_namelist_dir=$user_namelist_dir
init_data_dir=$init_data_dir

walltime="${walltime}"

export PROJECT=${project_code}

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
EOF

cat $wk_dir/env_settings.sh >> $script_file

cat << EOF >> $script_file

$cesm_create_newcase \
-case      \$casename                \
-compset   \$compset                     \
-res       \$resolution                       \
-mach      \$machine

cd \$casename

setXML "env_run.xml" "\${env_run[@]}"
setXML "env_mach_pes.xml" "\${env_mach_pes[@]}"

# Must setup here to get calculated TOTALPES
./cesm_setup

getXML "\${env_vars[@]}"

nodes=\$(( \$totalpes / \$max_tasks_per_node ))

# copy user namelist
cp \$user_namelist_dir/user_nl_* .



cat << XEOFX > config.jl

let
global overwrite_configs = Dict(
    "casename"                   => "\${casename}",
    "caseroot"                   => "\${caseroot}",
    "caserun"                    => "\${caserun}",
    "domain_file"                => "\${ocn_domain_path}/\${ocn_domain_file}",
    "short_term_archive_dir"     => "\${caserun}",
    "long_term_archive_dir"      => "\${dout_s_root}/ocn/hist",
    "enable_short_term_archive"  => true,
    "enable_long_term_archive"   => true,
    "daily_record"               => false,
    "monthly_record"             => true,
    "yearly_snapshot"            => true,
    "short_term_archive_list"    => "SSM_short_term_archive_list.txt",
    "substeps"                   => 8,
    "init_file"                  => "\${init_file}",
)
end

XEOFX

cat << XEOFX >> config.jl
$( cat $wk_dir/config_specific/config_${model}.jl )
XEOFX


cat << XEOFX > run_ocn.sh

#!/bin/bash
#PBS -N \$casename-ocn
#PBS -l walltime=\$walltime
#PBS -q share
### Merge output and error files
#PBS -j oe
#PBS -l select=1:ncpus=36:mpiprocs=36
### Send email on abort, begin and end
#PBS -m abe
#PBS -M meteorologytoday@gmail.com

### Run the ocean model ###

LID="\$(date +%y%m%d-%H%M)"
ocn_code="\$caseroot/SMARTSLAB-main/src/CESM_driver/run.jl"
config_file="\$caseroot/config.jl"

julia \$ocn_code --config="\$config_file" --core=${model} | tee -a SMARTSLAB.log.\$LID

XEOFX



mv \$casename.run \$casename.cesm.run

cat << XEOFX > \$casename.run

#!/bin/bash

qsub -A \$PROJECT -l walltime="\$walltime" \$caseroot/\$casename.ocn.run
qsub -A \$PROJECT -l walltime="\$walltime" \$caseroot/\$casename.cesm.run 

XEOFX

chmod +x \$casename.ocn.run
chmod +x \$casename.run

# Insert code
git clone https://github.com/meteorologytoday/SMARTSLAB-main.git
cd ./SourceMods/src.docn
ln -s ../../SMARTSLAB-main/src/CESM_driver/cesm1_docn_comp_mod.F90 ./docn_comp_mod.F90
ln -s ../../SMARTSLAB-main/src/CESM_driver/ProgramTunnel .




EOF

