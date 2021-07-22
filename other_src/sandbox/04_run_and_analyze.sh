#!/bin/bash

source 00_path.sh

./01_run.sh ${1} ${2}
./02_analysis.sh ${1}
