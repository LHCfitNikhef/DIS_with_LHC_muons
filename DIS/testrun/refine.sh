#!/bin/bash

for prefix in 'pwgPOWHEG+PYTHIA8-output'
do
    mergedata 4 $prefix-??.top -o $prefix-max.top
    mergedata 5 $prefix-??.top -o $prefix-min.top
done

for prefix in 'pwgPOWHEG+PYTHIA8-output'
do
    for type in '11' 'min' 'max' '12' '1H' '21' '22' 'H1' 'HH'
    do
      	for index in `grep index $prefix-$type.top | awk '{print $2}'`  
        do
            pastegnudata "[$index]" $prefix-$type.top > $prefix-$type-$index.dat
            sed -i 's/D/E/g' $prefix-$type-$index.dat
        done        
    done
done
