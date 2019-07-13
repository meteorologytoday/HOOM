#!/bin/bash

export wk_dir=$( dirname $0 )
#echo "wk_dir: $wk_dir"

lopts=(
    code-output-dir
    init-files-dir
    label
    resolution
    walltime
    data-clim-T-file
    data-clim-S-file
    topo-file
    varname-T
    varname-S
    varname-depth
    old-domain-file
    new-domain-file
    T-unit
    cesm-create-newcase
    user-namelist-dir
    compset
    machine
    project-code
    model
    model-config
)

source $wk_dir/getopt_helper.sh

cat << EOF

Users should be aware that domain file your are providing should
coincide with the ocean domain that is being used by CESM. This
is an inevitable pivot that cannot be overcome for now.

EOF

echo "Making directories..."
mkdir -p $code_output_dir
mkdir -p $init_files_dir


echo "Making initial files..."


new_data_clim_T_file=$init_files_dir/${label}_$( basename $data_clim_T_file ".nc" ).nc
new_data_clim_S_file=$init_files_dir/${label}_$( basename $data_clim_S_file ".nc" ).nc
new_topo_file=$init_files_dir/${label}_$( basename $topo_file ".nc" ).nc

zdomain_file=$init_files_dir/zdomain.nc

$wk_dir/make_init.sh                            \
    --output-dir=$init_files_dir                \
    --label=$label                              \
    --input-clim-T-file=$data_clim_T_file       \
    --input-clim-S-file=$data_clim_S_file       \
    --input-topo-file=$topo_file                \
    --output-clim-T-file=$new_data_clim_T_file  \
    --output-clim-S-file=$new_data_clim_S_file  \
    --output-topo-file=$new_topo_file           \
    --old-domain-file=$old_domain_file          \
    --new-domain-file=$new_domain_file          \
    --T-unit=$T_unit                            \
    --output-zdomain-file=$zdomain_file


echo "Making initial files for a specific model"
$wk_dir/make_init_each_model.sh                 \
    --output-dir=$init_files_dir                \
    --label=$label                              \
    --data-clim-T-file=$new_data_clim_T_file    \
    --data-clim-S-file=$new_data_clim_S_file    \
    --domain-file=$new_domain_file              \
    --zdomain-file=$zdomain_file                \
    --topo-file=$new_topo_file                      \
    --T-unit=$T_unit                            \
    --model=ESOM                                \
    --model-config=default


echo "Generate cesm sugar scripts..."
$wk_dir/make_cesm_sugar_script.sh           \
    --code-output-dir=$code_output_dir      \
    --resolution=$resolution                \
    --label=$label                          \
    --walltime="$walltime"                  \
    --project-code="$project_code"          
    


echo "Done."
