#!/bin/bash

wk_dir=$( dirname $0 )
echo $( basename $0 )
echo "wk_dir: $wk_dir"

lopts=(
    casename
    code-output-dir
    init-file
    label
    resolution
    walltime
    project-code
    compset
    machine
    cesm-create-newcase
    cesm-env
    user-namelist-dir
    model
    init-config
    ocn-ncpu
    qflux-file
    ocn-branch
    single-job
    machine
)

source $wk_dir/getopt_helper.sh

if [ ! -z "$casename" ]; then
    export casename="${label}_${resolution}_${model}_${init_config}"
fi

script_file=$code_output_dir/makecase_$casename.sh

# Clean the file
echo -ne "" > $script_file 

cat $wk_dir/lib_XML.sh >> $script_file

cat << EOF >> $script_file

casename=$casename
resolution=$resolution
machine=$machine
compset=$compset
user_namelist_dir=$user_namelist_dir
init_file=$init_file
qflux_file="$qflux_file"

walltime="${walltime}"
single_job="${single_job}"
ocn_ncpu=$ocn_ncpu

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

cat $cesm_env >> $script_file

cat << EOF >> $script_file

if [ -d \$casename ]; then
    echo "Error: \$casename already exists. Abort."
    exit 1;
fi



$cesm_create_newcase         \\
    -case      \$casename    \\
    -compset   \$compset     \\
    -res       \$resolution  \\
    -mach      \$machine

if [ ! -d \$casename ] ; then
    echo "Error: \$casename is not properly created. Abort."
    exit 1;
fi


cd \$casename


setXML "env_run.xml" "\${env_run[@]}"
setXML "env_mach_pes.xml" "\${env_mach_pes[@]}"

if [ -n "\$qflux_file" ]; then

    echo "Qflux file nonempty. Now setting user-defined qflux."
    setXML "env_run.xml" "DOCN_SOM_FILENAME" "\$qflux_file"
  
fi


# Must setup here to get calculated TOTALPES
./cesm_setup

getXML "\${env_vars[@]}"

nodes=\$(( \$totalpes / \$max_tasks_per_node ))

# copy user namelist
if [ "\$user_namelist_dir" != "" ]; then
    cp \$user_namelist_dir/user_nl_* .
fi

if [ ! -z "\$qflux_file" ]; then

    FORCING_DIR=\$( dirname \$qflux_file )
    FORCING_FILENAME=\$( basename \$qflux_file )

    cat << XEOFX > user_docn.streams.txt.som
    $( echo "$( cat $wk_dir/docn_stream.txt )" )
XEOFX

fi

cat << XEOFX > config.jl

let
global overwrite_configs = Dict(
    :casename                   => "\${casename}",
    :caseroot                   => "\${caseroot}",
    :caserun                    => "\${caserun}",
    :domain_file                => "\${ocn_domain_path}/\${ocn_domain_file}",
    :short_term_archive_dir     => "\${caserun}",
    :long_term_archive_dir      => "\${dout_s_root}/ocn/hist",
    :enable_short_term_archive  => true,
    :enable_long_term_archive   => true,
    :daily_record               => false,
    :monthly_record             => true,
    :yearly_snapshot            => true,
    :short_term_archive_list    => "SSM_short_term_archive_list.txt",
    :substeps                   => 8,
    :init_file                  => "\${init_file}",
)
end

XEOFX

cat << XEOFX >> config.jl
#  $wk_dir/config_specific/config_${model}.jl
$( cat $wk_dir/init_code/${model}_${init_config}/config.jl )
XEOFX


cat << XEOFX > \$casename.ocn.run

#!/bin/bash

### Run the ocean model ###

LID="\\\$(date +%y%m%d-%H%M)"
ocn_code="\$caseroot/SMARTSLAB-main/src/CESM_driver/run.jl"
config_file="\$caseroot/config.jl"
ocn_ncpu=\$ocn_ncpu

julia -p \\\$ocn_ncpu \\\$ocn_code --config="\\\$config_file" --core=${model} | tee -a SMARTSLAB.log.\\\$LID

XEOFX



mv \${casename}.run \${casename}.cesm.run

# ===== JOB SUBMISSION BLOCK BEGIN =====

if [ "\$single_job" != "on" ]; then

    echo "single_job != on. Use 2 jobs to complete a run."

    # "BATCHSUBMIT" is used when automatically resubmitting the job
    # Currently the design is to submit batch job through another script
    # file. So BATCHSUBMIT is set to just bash.

    setXML "env_run.xml" "BATCHSUBMIT" "bash"

    cat << XEOFX > \${casename}.run
$( echo "$( cat $wk_dir/batchsubmit_directives/two_jobs.${machine} )" )
XEOFX

else
    echo "single_job = on. Use 1 job to complete a run."

    if [ "\$nodes" -gt "1" ]; then
        echo "ERROR: Only 1 node is allowed when single_job variable is set as 'on'."
        exit 1
    fi

    # Single job Run (Experimental. Currently may only works on a single node because CESM 
    # was designed to take all the resources of each nodes. This scripts needs "env_mach_pes.xml"
    # configured correctly. If a single node got N cores. It would be M cores for cesm and (N-M)
    # cores for ocn model.
    
    cat << XEOFX > \${casename}.run
$( echo "$( cat $wk_dir/batchsubmit_directives/one_job.${machine} )" )
XEOFX
fi

# ===== JOB SUBMISSION BLOCK END =====

chmod +x \$casename.ocn.run
chmod +x \$casename.run

# Insert code
git clone --branch "$ocn_branch" https://github.com/meteorologytoday/SMARTSLAB-main.git

cd ./SourceMods/src.docn
ln -s ../../SMARTSLAB-main/src/CESM_driver/cesm1_docn_comp_mod.F90 ./docn_comp_mod.F90
ln -s ../../SMARTSLAB-main/src/CESM_driver/ProgramTunnel .


EOF

chmod +x $script_file 

echo "$casename"
