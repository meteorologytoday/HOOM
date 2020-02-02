#!/bin/bash



for m in $( seq 1 12 ); do

    echo "Doing month: $m"

    convert Q{1,2}_Qflux-$( printf "%02d" $m ).png +append tmp12.png
    convert Q{3,4}_Qflux-$( printf "%02d" $m ).png +append tmp34.png
    convert tmp12.png tmp34.png -append $( printf "COMBINE-Qflux-%02d.png" $m )

done





