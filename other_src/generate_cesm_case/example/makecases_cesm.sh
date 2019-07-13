#!/bin/bash

pwd_dir=$(pwd)

echo "Current directory: $pwd_dir"
machine=xtt-centos-intel
project_code=UCIR0029

T_file=$pwd_dir/raw_init_cond/b.e11.B1850C5CN.f09_g16.005.pop.h.TEMP.100001-109912.nc
S_file=$pwd_dir/raw_init_cond/b.e11.B1850C5CN.f09_g16.005.pop.h.SALT.100001-109912.nc

topo_file=$pwd_dir/CESM_domains/ocean_topog_gx1v6.nc

old_domain=$pwd_dir/CESM_domains/domain.ocn.gx1v6.090206.nc
new_domain=$pwd_dir/CESM_domains/domain.ocn.gx3v7.120323.nc

code_output_dir=$pwd_dir/cesm_scripts
init_files_dir=$pwd_dir/init_cond
cesm_env=$pwd_dir/env_settings.sh


model_settings=(
    SOM   default
    ESOM  default 
    MLMML restricted_climate  
    MLMML unrestricted       
)
#model_settings=(
#    ESOM  default
#)


for i in $(seq 1 $((${#model_settings[@]}/2))); do
    model=${model_settings[$((2*(i-1)))]}
    init_config=${model_settings[$((2*(i-1)+1))]}

    printf "Gonna making %s-%s\n" $model $model_config

    SMARTSLAB-code/other_src/generate_cesm_case/main.sh \
        --code-output-dir=$code_output_dir              \
        --init-files-dir=$init_files_dir                \
        --label=LENS                                    \
        --resolution=f45_g37                            \
        --walltime="06:00:00"                           \
        --data-clim-T-file=$T_file                      \
        --data-clim-S-file=$S_file                      \
        --topo-file=$topo_file                          \
        --old-domain-file=$old_domain                   \
        --new-domain-file=$new_domain                   \
        --T-unit=C                                      \
        --cesm-create-newcase=~/ucar_models/cesm1_2_2_1/scripts/create_newcase \
        --compset=E_1850_CN                             \
        --machine=$machine                              \
        --model=$model                                  \
        --init-config=$init_config                      \
        --cesm-env=$cesm_env                            \
        --project-code=$project_code


done
