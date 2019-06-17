#!/bin/bash

if [ -f $ocn_outputfile ]; then

    echo "No need to concat ocn files."
    echo "$ocn_outputfile already exists."

else

    printf "Concat ocn files... "
        # atm variables
    cd $ocn_hist_path 
    eval $(cat <<EOF
    ncrcat -O $casename.ocn.h.monthly.{$beg_year..$end_year}.nc $ocn_outputfile

EOF
    )
    cd $wpath

    printf "done.\n"

fi
