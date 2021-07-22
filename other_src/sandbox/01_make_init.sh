#!/bin/bash

source 00_path.sh

julia SMARTSLAB-main/other_src/CoordTrans/mk_HOOM_zdomain.jl --output-file="$zdomain_file" --resolution=Standard 
julia make_init.jl --output-file="$init_file" --domain-file="$domain_file" --topo-file="$topo_file" --zdomain-file="$zdomain_file"

