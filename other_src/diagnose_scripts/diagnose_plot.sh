#!/bin/bash

plot_path=$diagnose_scripts_path/plot

# atm
python3 $plot_path/plot_SST.py --data-file=$atm_outputfile --domain-file=$atm_domain --output-dir=$diagnose_output_dir

# ocn
python3 $plot_path/plot_ocean_diagnose.py --data-file-SSTAYYC=$ocn_trans_outputfile_mstat --data-file-SSTAVAR=$ocn_trans_outputfile_mstat --domain-file=$atm_domain --output-dir=$diagnose_output_dir




# General
python3 $plot_path/plot_climate_indices.py --data-file-PDO=$ocn_trans_outputfile_SSTA --output-dir=$diagnose_output_dir --data-file-AO=$atm_outputfile_anomalies

